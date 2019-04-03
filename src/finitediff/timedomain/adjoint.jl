"""
   The adjoint operator of one_step_forward.
"""
function one_step_adjoint!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil)

    spt2.vz = spt1.vz + (fidMtx.MpzBvz)' * spt1.pz
    spt2.vx = spt1.vx + (fidMtx.MpxBvx)' * spt1.px
    spt2.pz =            fidMtx.MpzBpz  .* spt1.pz
    spt2.px =            fidMtx.MpxBpx  .* spt1.px

    spt2.pz = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.pz
    spt2.px = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.px
    spt2.vz =  fidMtx.MvzBvz  .* spt2.vz
    spt2.vx =  fidMtx.MvxBvx  .* spt2.vx

    return nothing
end

"""
   The adjoint operator of one_step_forward (in-place).
"""
function one_step_adjoint!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil,
         tmp1::Vector{Tv}, tmp2::Vector{Tv}) where {Tv<:AbstractFloat}

    At_mul_b!(tmp1, ofds.MpzBvz, spt1.pz)
    addition!(spt2.vz, spt1.vz, tmp1)

    At_mul_b!(tmp1, ofds.MpxBvx, spt1.px)
    addition!(spt2.vx, spt1.vx, tmp1)

    multiplication!(spt2.pz, ofds.MpzBpz, spt1.pz)
    multiplication!(spt2.px, ofds.MpxBpx, spt1.px)

    At_mul_b!(tmp1, ofds.MvzBp, spt2.vz)
    At_mul_b!(tmp2, ofds.MvxBp, spt2.vx)
    addition!(tmp1, tmp2)

    addition!(spt2.pz, tmp1)
    addition!(spt2.px, tmp1)

    multiplication!(spt2.vz, ofds.MvzBvz)
    multiplication!(spt2.vx, ofds.MvxBvx)

    return nothing
end

"""
   one step backward reconstruction of wave field from the boundary values
"""
function one_step_backward!(wfd1::Wavefield, wfd2::Wavefield, rfds::RigidFDStencil, it::Ti,
         bnd::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, params::ModelParams) where {Ti<:Int64, Tv<:AbstractFloat}

    # wfd1.p  = wfd2.p  - rfds.MpxBvx * wfd2.vx - rfds.MpzBvz * wfd2.vz
    # wfd1.p[bound.index]  = bound.p[:, it]
    A_mul_b!(tmp1, rfds.MpxBvx, wfd2.vx)
    A_mul_b!(tmp2, rfds.MpzBvz, wfd2.vz)
    addition!(tmp1, tmp2)
    minus!(wfd1.p, wfd2.p, tmp1)

    # correct for boundary part of pressure
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd1.p[j] = bnd.p[i,it]
    end

    # wfd1.vx = wfd2.vx - rfds.MvxBp  * wfd1.p
    # wfd1.vz = wfd2.vz - rfds.MvzBp  * wfd1.p
    # wfd1.vx[bound.index] = bound.vx[:, it]
    # wfd1.vz[bound.index] = bound.vz[:, it]
    A_mul_b!(tmp1, rfds.MvxBp, wfd1.p)
    A_mul_b!(tmp2, rfds.MvzBp, wfd1.p)
    minus!(wfd1.vx, wfd2.vx, tmp1)
    minus!(wfd1.vz, wfd2.vz, tmp2)

    # correct for boundary part of particle velocity
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd1.vz[j] = bnd.vz[i,it]
        wfd1.vx[j] = bnd.vx[i,it]
    end

    return nothing
end

"""
   the adjoint operator of one point source simulation.
"""
function multi_step_adjoint(rec::Recordings, ofds::ObsorbFDStencil,
         src::Source, params::ModelParams)

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # inject the recordings to the last snapshot
    inject_rec2spt!(spt2, rec, params.nt)

    # the time range of source function
    p   = zeros(params.data_format, src.it_max - src.it_min + 1)
    idx = src.src2spt

    if src.it_min <= params.nt <= src.it_max
       i = params.nt - src.it_min + 1
       p[i] = spt2.pz[idx] * params.dt
    end

    # back propagation
    for it = params.nt-1 : -1 : 1

        one_step_adjoint!(spt1, spt2, ofds, tmp1, tmp2)
        inject_rec2spt!(spt1, rec, it)

        if src.it_min <= it <= src.it_max
           i = it - src.it_min + 1
           p[i] = spt1.pz[idx] * params.dt
        end

        copy_snapshot!(spt2, spt1)
    end

    return p
end

"""
   Reconstruct pressure field backward using the boundary wavefield value and
the last wavefield
"""
function pressure_reconstruct_backward(bnd::WavefieldBound, wfd::Wavefield,
         rfds::RigidFDStencil, src::Source, params::ModelParams)

    # length of one-step pressure field
    N = params.nz * params.nx

    # decide the number of boundary layers
    order = length(params.fd_coefficients)

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = Wavefield(params)
    tmp1 = zeros(params.data_format, N)
    tmp2 = zeros(params.data_format, N)

    copy_wavefield!(wfd2, wfd)

    # pre-allocate memory for saving the reconstructed pressure field
    pre = zeros(params.data_format, N * params.nt)
    idx_o = N*params.nt - N + 1 # the pressure field at last time step

    for it = params.nt : -1 : 1

        # one step back
        one_step_backward!(wfd1, wfd2, rfds, it, bnd, tmp1, tmp2, params)

        # save the pressure field
        copyto!(pre, idx_o, wfd1.p, 1, N)
        idx_o = idx_o - N

        # prepare for next step
        copy_wavefield!(wfd2, wfd1)
        subtract_source!(wfd2, src, it)
    end

    return reshape(pre, params.nz, params.nx, params.nt)
end
