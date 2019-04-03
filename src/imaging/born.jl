"""
   compute the time derivative of pressure, which is savd as a matrix, the second
axis is time step. The time derivative is approximated by central finite-difference.
"""
function compute_dpdt!(dpdt::Matrix{Tv}, spt1::Snapshot, spt2::Snapshot,
         mp::Vector{Tv}, spt2wfd::Vector{Ti}, it::Ti, double_dt::Tv) where {Ti<:Int64, Tv<:AbstractFloat}

    for i = 1 : length(spt2wfd)
        # mapping spt to pressure
        j = spt2wfd[i]

        # central finite difference
        dpdt[i,it-1] = (spt2.pz[j]+spt2.px[j]-mp[i]) / double_dt

        # update memory variable
        mp[i] = spt1.pz[j] + spt1.px[j]
    end
    return nothing
end

"""
   compute the time derivative of source-side wavefield at current time step via central
finite difference.
"""
function compute_one_dpdt!(dpdt::Vector{Tv}, spt1::Snapshot, spt2::Snapshot,
         mp::Vector{Tv}, params::ModelParams) where {Tv<:AbstractFloat}

    N = params.nz * params.nx
    one_over_2dt = params.data_format(1.0 / (2.0 * params.dt))

    for i = 1 : N

        # mapping spt to pressure
        j = params.spt2wfd[i]

        # central finite difference
        dpdt[i] = (spt2.pz[j]+spt2.px[j]-mp[i]) * one_over_2dt

        # update memory variable
        mp[i] = spt1.pz[j] + spt1.px[j]
    end
    return nothing
end

"""
   compute the source-side wavefield dpdt, time derivative is approximated by central finite-difference.
"""
function get_sourceside_wavefield(src::Source, ofds::ObsorbFDStencil, params::ModelParams)

    # initialize some variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # 2 * dt, used for time derivative approximation
    double_dt = params.data_format(2.0 * params.dt)

    # initialize the a matrix to save the result
    N = params.nz * params.nx
    dpdt = zeros(params.data_format, N, params.nt)

    # memory variable save the pressure field at it-1
    mp   = zeros(params.data_format, N)

    # add source
    add_source!(spt1, src, 1)

    # time stepping
    for it = 2 : params.nt+1    #to compute dpdt at nt step, we need the pressure field at nt-1 and nt+1

        # one time step forward
        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)

        # compute virtual source
        compute_dpdt!(dpdt, spt1, spt2, mp, params.spt2wfd, it, double_dt)

        # prepare for next step
        copy_snapshot!(spt1, spt2)
    end

    # reshape into a cube
    return reshape(dpdt, params.nz, params.nx, params.nt)
end

"""
   Add the virtual source to pressure (just add to pz), the source-side wavefield
of all time-steps are saved in memory.
"""
function add_virtual_source!(spt::Snapshot, dpdt::Matrix{Tv}, delta_lambda::Vector{Tv},
         params::ModelParams, it::Ti) where {Tv<:AbstractFloat, Ti<:Int64}

    N = params.nz * params.nx
    for i = 1 : N
        j = params.spt2wfd[i]
        spt.pz[j] = spt.pz[j] + dpdt[i,it] * params.dt * delta_lambda[i]
    end
    return nothing
end

"""
   Add the virtual source to pressure field, the source-side wavefield at current time step is
computed on the flight
"""
function add_virtual_source!(spt::Snapshot, dpdt::Vector{Tv}, delta_lambda::Vector{Tv},
         params::ModelParams) where {Tv<:AbstractFloat}

    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]
        spt.pz[j] = spt.pz[j] + dpdt[i] * params.dt * delta_lambda[i]
    end
    return nothing
end

"""
   forward born approximation with source side wavefield saved in memory
"""
function born_approximation_forward!(rec::Recordings, dpdt::Matrix{Tv}, delta_lambda::Vector{Tv},
         ofds::ObsorbFDStencil, params::ModelParams) where{Ti<:Int64, Tv<:AbstractFloat}

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # add virtual source to the first snapshot
    add_virtual_source!(spt1, dpdt, delta_lambda, params, 1)
    sample_spt2rec!(rec, spt1, 1)

    # loop over time step
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_virtual_source!(spt2, dpdt, delta_lambda, params, it)
        sample_spt2rec!(rec, spt2, it)

        # prepare for next time step
        copy_snapshot!(spt1, spt2)
    end

    return rec
end

"""
   forward Born approximation with source-side wavefield reconstruction are computed on the flight
"""
function born_approximation_forward!(rec::Recordings, delta_lambda::Vector{Tv}, src::Source,
         ofds::ObsorbFDStencil, params::ModelParams) where{Tv<:AbstractFloat}

    # length of model
    N = params.nz * params.nx

    # intermediate variables born forward time stepping
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # intermediate variables for computing source side wavefield
    spt_s1 = Snapshot(params)
    spt_s2 = Snapshot(params)

    # memory variable to save the pressure field at previous time step
    dpdt = zeros(params.data_format, N)
    mp   = zeros(params.data_format, N)
    cube = zeros(params.data_format, N, params.nt)
    idx_o= 1

    # add source to wavefield to get the first pressure field
    add_source!(spt_s1, src, 1)

    for it = 1 : params.nt

        # compute the source side wavefield (it+1)
        one_step_forward!(spt_s2, spt_s1, ofds, tmp1, tmp2)
        add_source!(spt_s2, src, it+1)

        # compute dpdt at it
        compute_one_dpdt!(dpdt, spt_s1, spt_s2, mp, params)
        copyto!(cube, idx_o, dpdt, 1, N);
        idx_o = idx_o+N

        # add virtual source
        add_virtual_source!(spt1, dpdt, delta_lambda, params)

        # sample the scattered wave field at it
        sample_spt2rec!(rec, spt1, it)

        # prepare for source-side wavefield at (it+2)
        copy_snapshot!(spt_s2, spt_s1)

        # the scattered wavefield at it+1 without adding virtual source
        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)

        # prepare for adding virtual source at it+1
        copy_snapshot!(spt1, spt2)
    end

    return cube
end

"""
   apply imaging condition to collapse receiver-side wavefield, which is the adjoint
operator of add virtual sources
"""
function apply_image_condition!(delta_lambda::Vector{Tv}, spt::Snapshot,
         dpdt::Matrix{Tv}, params::ModelParams, it::Ti) where {Ti<:Int64, Tv<:AbstractFloat}

    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]
        delta_lambda[i] = delta_lambda[i] + spt.pz[j] * dpdt[i,it] * params.dt
    end
    return nothing
end

"""
   the adjoint operator of add virtual source, where the source-side wavefield is computed on the flight
"""
function apply_image_condition!(delta_lambda::Vector{Tv}, spt::Snapshot, dpdt::Vector{Tv},
         params::ModelParams) where {Tv<:AbstractFloat}

    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]
        delta_lambda[i] = delta_lambda[i] + spt.pz[j] * dpdt[i] * params.dt
    end
    return nothing
end

"""
   The Adjoint operator of Born approximation with source side wavefield saved in memory
"""
function born_approximation_adjoint(rec::Recordings, dpdt::Matrix{Tv},
         ofds::ObsorbFDStencil, params::ModelParams) where {Tv<:AbstractFloat}

    # length of gradient
    N = params.nz * params.nx

    # pre-allocate memory for gradient
    delta_lambda = zeros(params.data_format, N)

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # inject recordings to the last snapshot
    inject_rec2spt!(spt2, rec, params.nt)
    apply_image_condition!(delta_lambda, spt2, dpdt, params, params.nt)

    for it = params.nt-1 : -1 : 1

        one_step_adjoint!(spt1, spt2, ofds, tmp1, tmp2)
        inject_rec2spt!(spt1, rec, it)

        apply_image_condition!(delta_lambda, spt1, dpdt, params, it)

        copy_snapshot!(spt2, spt1)
    end

    return delta_lambda
end

# function born_approximation_adjoint(rec::Record, fidMtx::FiniteDiffMatrix, src::Source,
#          bnd::WavefieldBound, wfd::Wavefield, fidMtxT::RigidFiniteDiffMatrix, par::PhysicalModel)
#
#     N = par.nz * par.nx
#     const one_over_two = par.data_format(0.5)
#     delta_lambda = zeros(par.data_format, N)
#
#     spt1 = Snapshot(par)
#     spt2 = Snapshot(par)
#     tmp1 = zeros(spt1.vz)
#     tmp2 = zeros(spt1.vx)
#
#     wfd1 = Wavefield(par)
#     wfd2 = Wavefield(par)
#     tmp3 = zeros(wfd1.p)
#     tmp4 = zeros(wfd1.p)
#     mp   = copy(wfd.p)
#
#     inject_record2spt!(spt2, rec, par.nt)
#     one_step_backward!(wfd2, wfd , fidMtxT, par.nt  , bnd, tmp3, tmp4, par)
#     one_step_backward!(wfd1, wfd2, fidMtxT, par.nt-1, bnd, tmp3, tmp4, par)
#     get_time_derivative!(tmp3, mp, wfd1.p, one_over_two)
#     apply_image_condition!(delta_lambda, spt2, tmp3, par)
#
#     for it = par.nt-1 : -1 : 2
#
#         # backward propagation of recordings
#         one_step_adjoint!(spt1, spt2, fidMtx, tmp1, tmp2)
#         inject_record2spt!(spt1, rec, it)
#
#         # update memory variable
#         Base.LinAlg.BLAS.blascopy!(N, wfd2.p, 1, mp, 1)
#         if it <= src.nt-1
#            mp[src.idx2] = mp[src.idx2] + src.p[it+1]*par.dt
#         end
#
#         # reconstruct of source-side wavefield
#         copy_wavefield!(wfd2, wfd1)
#         if it <= src.nt
#            wfd2.p[src.idx2] = wfd2.p[src.idx2] - src.p[it]*par.dt
#         end
#         one_step_backward!(wfd1, wfd2, fidMtxT, it-1, bnd, tmp3, tmp4, par)
#
#         # apply imaging condition
#         get_time_derivative!(tmp3, mp, wfd1.p, one_over_two)
#         apply_image_condition!(delta_lambda, spt1, tmp3, par)
#
#         copy_snapshot!(spt2, spt1)
#     end
#
#     # the last one
#     one_step_adjoint!(spt1, spt2, fidMtx, tmp1, tmp2)
#     inject_record2spt!(spt1, rec, 1)
#
#     wfd2.p[src.idx2] = wfd2.p[src.idx2] + src.p[2] * par.dt
#     for i = 1 : N
#         j = par.index1[i]
#         delta_lambda[i] = delta_lambda[i] + spt1.pz[j] * wfd2.p[i] * one_over_two
#     end
#
#     return delta_lambda
# end
#
# # one shot reverse time migration
# function single_shot_RTM(rec::Record, fidMtx::FiniteDiffMatrix, src::Source,
#          bnd::WavefieldBound, wfd::Wavefield, fidMtxT::RigidFiniteDiffMatrix, par::PhysicalModel)
#
#     N = par.nz * par.nx
#     const one_over_two = par.data_format(0.5)
#     delta_lambda = zeros(par.data_format, N)
#     src_side     = zeros(par.data_format, N)
#     rec_side     = zeros(par.data_format, N)
#
#     spt1 = Snapshot(par)
#     spt2 = Snapshot(par)
#     tmp1 = zeros(spt1.vz)
#     tmp2 = zeros(spt1.vx)
#
#     wfd1 = Wavefield(par)
#     wfd2 = Wavefield(par)
#     tmp3 = zeros(wfd1.p)
#     tmp4 = zeros(wfd1.p)
#     mp   = copy(wfd.p)
#
#     inject_record2spt!(spt2, rec, par.nt)
#     one_step_backward!(wfd2, wfd , fidMtxT, par.nt  , bnd, tmp3, tmp4, par)
#     one_step_backward!(wfd1, wfd2, fidMtxT, par.nt-1, bnd, tmp3, tmp4, par)
#     get_time_derivative!(tmp3, mp, wfd1.p, one_over_two)
#     apply_normalized_image_condition!(delta_lambda, src_side, rec_side, spt2, tmp3, par)
#
#     for it = par.nt-1 : -1 : 2
#
#         # backward propagation of recordings
#         one_step_adjoint!(spt1, spt2, fidMtx, tmp1, tmp2)
#         inject_record2spt!(spt1, rec, it)
#
#         # update memory variable
#         Base.LinAlg.BLAS.blascopy!(N, wfd2.p, 1, mp, 1)
#         if it <= src.nt-1
#            mp[src.idx2] = mp[src.idx2] + src.p[it+1]*par.dt
#         end
#
#         # reconstruct of source-side wavefield
#         copy_wavefield!(wfd2, wfd1)
#         if it <= src.nt
#            wfd2.p[src.idx2] = wfd2.p[src.idx2] - src.p[it]*par.dt
#         end
#         one_step_backward!(wfd1, wfd2, fidMtxT, it-1, bnd, tmp3, tmp4, par)
#
#         # apply imaging condition
#         get_time_derivative!(tmp3, mp, wfd1.p, one_over_two)
#         apply_normalized_image_condition!(delta_lambda, src_side, rec_side, spt1, tmp3, par)
#
#         copy_snapshot!(spt2, spt1)
#     end
#
#     # the last one
#     one_step_adjoint!(spt1, spt2, fidMtx, tmp1, tmp2)
#     inject_record2spt!(spt1, rec, 1)
#
#     wfd2.p[src.idx2] = wfd2.p[src.idx2] + src.p[2] * par.dt
#     for i = 1 : N
#         j = par.index1[i]
#         delta_lambda[i] = delta_lambda[i] + spt1.pz[j] * wfd2.p[i] * one_over_two
#         src_side[i]     = src_side[i] + (wfd2.p[i] * one_over_two)^2
#         rec_side[i]     = rec_side[i] + (spt1.pz[i])^2
#     end
#
#     return delta_lambda, src_side, rec_side
# end
#
# """
#    correct the wavefield by source
# """
# function correct_wavefield!(wfd::Wavefield, src::Source, it::Int64, par::PhysicalModel)
#
#     tc = (it-1) * src.dt
#     if src.ot <= tc <= src.tmax && par.order+1 <= src.isz <= par.nz-par.order && par.order+1 <= src.isx <= par.nx-par.order
#        indt = it - round(Int64, src.ot / src.dt)
#        pos  = src.idx2
#        wfd.p[pos] = wfd.p[pos] + src.p[indt] * src.dt
#     end
#     return nothing
# end


# """
#    compute the frequency component of source-side wavefield on the fly
# """
# function compute_dpdt_dft!(dpdt::Array{Complex{Tv},3}, spt1::Snapshot, spt2::Snapshot,
#          mp::Vector{Tv}, basis::Vector{Complex{Tv}}, factor::Vector{Complex{Tv}}, nf::Ti, one_over_2dt::Tv,
#          phi::PhysicalModel) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     if phi.data_format == Float32
#
#        ccall((:compute_dpdt_dft_float32_, spmatveclib), Void,
#              (Ref{Complex64}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Int64}, Ref{Complex64},
#               Ref{Complex64}, Ref{Float32}, Ref{Int64}  , Ref{Int64}  , Ref{Int64}) ,
#               dpdt          , spt1.pz     , spt1.px     , spt2.pz     , spt2.px     , mp          , phi.index1, basis         ,
#               factor        , one_over_2dt, nf          , phi.nz      , phi.nx       )
#
#     elseif phi.data_format == Float64
#
#        ccall((:compute_dpdt_dft_float64_, spmatveclib), Void,
#              (Ref{Complex128}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Complex128},
#               Ref{Complex128}, Ref{Float64}, Ref{Int64}  , Ref{Int64}  , Ref{Int64}) ,
#               dpdt           , spt1.pz     , spt1.px     , spt2.pz     , spt2.px     , mp          , phi.index1, basis          ,
#               factor         , one_over_2dt, nf          , phi.nz      , phi.nx       )
#
#     end
#     return nothing
# end
