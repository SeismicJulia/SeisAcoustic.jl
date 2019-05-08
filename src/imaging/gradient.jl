# """
#    compute the source-side wave field as (p[it+1/2]-p[it-1/2]-f[it]) / k
# """
# function apply_image_condition!(g::Vector{Tv}, p_adj::Vector{Tv}, wfd1::Wavefield,
#          wfd2::Wavefield, params::ModelParams) where {Tv <: AbstractFloat}
#
#     # number of elements
#     N = params.nz * params.nx
#
#     # constant
#     # double = params.data_format(2.0)
#
#     # compute gradient
#     for i = 1 : N
#         # dpdt = (wfd2.p[i] - wfd1.p[i]) / (params.rho[i] * params.vel[i] * params.vel[i])
#         dpdt = (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
#         g[i] = g[i] + dpdt * p_adj[i]
#     end
#     return nothing
# end
#
# """
#    compute the gradient of bulk modulus
# """
# function bulk_gradient(rec::Recordings, path_bnd::Ts, path_wfd::Ts,
#          src::Source, params::ModelParams) where {Ts <: String}
#
#     # model length
#     N = params.nz * params.nx
#
#     # allocate memory for computing adjoint wavefield
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp  = zeros(params.data_format, params.Nz * params.Nx)
#     tmp_z1 = zeros(params.data_format, params.Nz)
#     tmp_z2 = zeros(params.data_format, params.Nz)
#     tmp_x1 = zeros(params.data_format, params.Nx)
#     tmp_x2 = zeros(params.data_format, params.Nx)
#
#     # allocate memory for reconstructing source-side wavefield backward
#     wfd1 = Wavefield(params)
#     wfd2 = read_one_wavefield(path_wfd, 1)
#     wfd_z1 = zeros(params.data_format, params.nz)
#     wfd_z2 = zeros(params.data_format, params.nz)
#     wfd_x1 = zeros(params.data_format, params.nx)
#     wfd_x2 = zeros(params.data_format, params.nx)
#
#     # ============================================
#     # pre = zeros(params.data_format, N * params.nt)
#     # idx_o = N*params.nt - N + 1
#     # copyto!(pre, idx_o, wfd2.p, 1, N)
#     # ============================================
#
#     # allocate memory to save the adjoint pressure field
#     p_adj = zeros(params.data_format, N)
#
#     # allocate memory to save the gradient
#     g     = zeros(params.data_format, N)
#
#     # inject the recordings to the last snapshot
#     inject_rec2spt!(spt2, rec, params.nt)
#
#     # get the adjoint wavefield at the last time step
#     sample_adjoint2pre!(p_adj, spt2, params)
#
#     # open the file of boundary wavefield
#     fid_bnd = open(path_bnd, "r")
#     bnd     = WavefieldBound(params)
#
#     # prepare for backward reconstruction
#     subtract_source!(wfd2, src, params.nt)
#     read_one_boundary!(bnd, fid_bnd, params.nt-1, params)
#     one_step_backward!(wfd1, wfd2, bnd, params,
#                        wfd_z1, wfd_z2, wfd_x1, wfd_x2)
#
#     # =================================
#     # idx_o = idx_o - N
#     # copyto!(pre, idx_o, wfd1.p, 1, N)
#     # =================================
#
#     # compute the gradient
#     apply_image_condition!(g, p_adj, wfd1, wfd2, params)
#
#     # prepare for the next step backward reconstruction
#     copy_wavefield!(wfd2, wfd1)
#
#     # backward propagation
#     for it = params.nt-1 : -1 : 2
#
#         # compute adjoint wavefield
#         one_step_adjoint!(spt1, spt2, params,
#                           tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
#         inject_rec2spt!(spt1, rec, it)
#         copy_snapshot!(spt2, spt1)
#
#         # get the adjoint pressure at the current time step
#         sample_adjoint2pre!(p_adj, spt2, params)
#
#         # read the source-side wavefield at the last time step
#         subtract_source!(wfd2, src, it)
#         read_one_boundary!(bnd, fid_bnd, it-1, params)
#         one_step_backward!(wfd1, wfd2, bnd, params,
#                            wfd_z1, wfd_z2, wfd_x1, wfd_x2)
#
#         # =================================
#         # idx_o = idx_o - N
#         # copyto!(pre, idx_o, wfd1.p, 1, N)
#         # =================================
#
#         # compute the gradient
#         apply_image_condition!(g, p_adj, wfd1, wfd2, params)
#
#         # prepare for the next step backward reconstruction
#         copy_wavefield!(wfd2, wfd1)
#     end
#
#     close(fid_bnd)
#     return g
# end

"""
   compute the source-side wave field by backward reconstruction with saved boundary values
"""
function source_side_wavefield(path_bnd::Ts, path_wfd::Ts,
         src::Source, params::ModelParams) where {Ts <: String}

    # model length
    N = params.nz * params.nx

    # allocate memory for reconstructing source-side wavefield backward
    wfd1 = Wavefield(params)
    wfd2 = read_one_wavefield(path_wfd, 1)
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # open the file of boundary wavefield
    fid_bnd = open(path_bnd, "r")
    bnd     = WavefieldBound(params)

    # prepare for backward reconstruction
    subtract_source!(wfd2, src, params.nt)
    read_one_boundary!(bnd, fid_bnd, params.nt-1, params)
    one_step_backward!(wfd1, wfd2, bnd, params,
                       wfd_z1, wfd_z2, wfd_x1, wfd_x2)

    # save nt-1 source side wavefield
    pre = zeros(params.data_format, N * (params.nt-1))
    tmp = zeros(params.data_format, N)
    idx_o = N * (params.nt-2) + 1

    for i = 1 : N
        tmp[i] = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
    end
    copyto!(pre, idx_o, tmp, 1, N)

    # prepare for the next step backward reconstruction
    copy_wavefield!(wfd2, wfd1)

    # backward propagation
    for it = params.nt-1 : -1 : 2

        # read the source-side wavefield at the last time step
        subtract_source!(wfd2, src, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)
        one_step_backward!(wfd1, wfd2, bnd, params,
                           wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        for i = 1 : N
            tmp[i] = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        end
        idx_o = idx_o - N
        copyto!(pre, idx_o, tmp, 1, N)

        # prepare for the next step backward reconstruction
        copy_wavefield!(wfd2, wfd1)
    end

    return reshape(pre, params.nz, params.nx, params.nt-1)
end

"""
   compute source-side wavefield by forward modeling, we provide two options to compute
the source-side wavefield. iflag=1: source side wavefield is computed from pressure
                           iflag=2: source side wavefield is computed from particle velocity
"""
function get_sourceside_wavefield(src::Source, params::ModelParams; iflag=1)

    N  = params.nz * params.nx
    tmp= zeros(params.data_format, N)

    # initialize variables for time stepping
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # save nt-1 source side wavefield
    pre = zeros(params.data_format, N * (params.nt-1))
    tmp = zeros(params.data_format, N)
    if iflag == 2
       rz = zeros(params.data_format, params.Nz * params.Nx)
       rx = zeros(params.data_format, params.Nz * params.Nx)
    end
    idx_o = 1

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # loop over time stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # save the source-side wavefield before adding source term to it
        if iflag == 1
           for i = 1 : N
               j = params.spt2wfd[i]
               tmp[i] = 2.0 * (spt2.pz[j] + spt2.px[j] - spt1.pz[j] - spt1.px[j]) / params.vel[i]
           end
        elseif iflag == 2
           # vertical partial derivative of vz
           vertical_partial_derivative!(rz, spt2.vz, params.dvdz, params.Nz, params.Nx,
                                        tmp_z1, tmp_z2)
           # horizontal partial derivative of vx
           horizontal_partial_derivative!(rx, spt2.vx, params.dvdx, params.Nz, params.Nx,
                                          tmp_x1, tmp_x2)
           # compute the source-side wavefield at one-time step
           for i = 1 : N
               j = params.spt2wfd[i]
               tmp[i] = 2.0 * params.rho[i] * params.vel[i] * params.dt * (rz[j]+rx[j])
           end
        end

        copyto!(pre, idx_o, tmp, 1, N)
        idx_o = idx_o + N

        # add the source term to it
        add_source!(spt2, src, it)

        # prepare for next iteration
        copy_snapshot!(spt1, spt2)
    end

    # reshape source-side wavefield into a cube
    pre =  reshape(pre, params.nz, params.nx, params.nt-1)

    # set the top to 0
    if iflag == 2
       pre[1,:,:] .=0.0
    end

    return pre
end
