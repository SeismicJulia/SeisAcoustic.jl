"""
   compute the residues between two recordings (synthetic and observed data)
"""
function get_residue(dsyn::Recordings, dobs::Recordings)

    # data format of recordings
    data_format = eltype(dsyn.p)

    # generate matrix to save the difference between dsyn and dobs
    p = zeros(data_format, dsyn.nt, dsyn.nr)

    # compute the difference
    for ir = 1 : dsyn.nr
        for it = 1 : dsyn.nt
            p[it,ir] = dsyn.p[it,ir] - dobs.p[it,ir]
        end
    end

    # return the residue in a format of recordings
    return Recordings(dsyn.nt, dsyn.nr, dsyn.dt, copy(dsyn.irz), copy(dsyn.irx), copy(dsyn.spt2rec), p)
end

"""
   compute the source-side wave field as (p[it+1/2]-p[it-1/2]-f[it]) / k
"""
function apply_image_condition!(g::Vector{Tv}, p_adj::Vector{Tv}, wfd1::Wavefield,
         wfd2::Wavefield, params::TdParams) where {Tv <: AbstractFloat}

    # number of elements
    N = params.nz * params.nx

    # compute gradient
    for i = 1 : N
        dpdt = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        g[i] = g[i] + dpdt * p_adj[i]
    end
    return nothing
end

"""
   compute the gradient of velocity model, the source-side wavefield is reconstructed from
last wavefield and boundary value, the adjoint wavefield is obtained by backward propagating
the residues.
"""
function velocity_gradient(res::Recordings, path_bnd::Ts, path_wfd::Ts,
         src::Source, params::TdParams) where {Ts <: String}

    # model length
    N = params.nz * params.nx

    # allocate memory for computing adjoint wavefield
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp    = zeros(params.data_format, params.Nz * params.Nx)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # allocate memory for reconstructing source-side wavefield backward
    wfd1 = Wavefield(params)
    wfd2 = read_one_wavefield(path_wfd, 1)
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # allocate memory to save the adjoint pressure field
    p_adj = zeros(params.data_format, N)

    # allocate memory to save the gradient
    g     = zeros(params.data_format, N)

    # inject the recordings to the last snapshot
    inject_rec2spt!(spt2, res, params.nt)

    # get the adjoint pressue wavefield at the last time step
    sample_adjoint2pre!(p_adj, spt2, params)

    # open the file of boundary wavefield
    fid_bnd = open(path_bnd, "r")
    bnd     = WavefieldBound(params)

    # backward reconstruct wavefield
    subtract_source!(wfd2, src, params.nt)
    read_one_boundary!(bnd, fid_bnd, params.nt-1, params)
    one_step_backward!(wfd1, wfd2, bnd, params,
                       wfd_z1, wfd_z2, wfd_x1, wfd_x2)

    # compute the gradient
    apply_image_condition!(g, p_adj, wfd1, wfd2, params)

    # prepare for the next step backward reconstruction
    copy_wavefield!(wfd2, wfd1)

    # backward propagation
    for it = params.nt-1 : -1 : 2

        # compute adjoint wavefield
        one_step_adjoint!(spt1, spt2, params,
                          tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, res, it)
        copy_snapshot!(spt2, spt1)

        # get the adjoint pressure at the current time step
        sample_adjoint2pre!(p_adj, spt2, params)

        # read the source-side wavefield at the last time step
        subtract_source!(wfd2, src, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)
        one_step_backward!(wfd1, wfd2, bnd, params,
                           wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # compute the gradient
        apply_image_condition!(g, p_adj, wfd1, wfd2, params)

        # prepare for the next step backward reconstruction
        copy_wavefield!(wfd2, wfd1)
    end

    # close the boundary value file
    close(fid_bnd)

    # return the gradient of velocity model
    return g
end

"""
   compute source-side wavefield by forward modeling, we provide two options to compute
the source-side wavefield. iflag=1: source side wavefield is computed from pressure
                           iflag=2: source side wavefield is computed from particle velocity
"""
function get_sourceside_wavefield(src::Source, params::TdParams; iflag=1)

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
    if params.free_surface && iflag == 2
       pre[1,:,:] .= 0.0
    end

    return pre
end
