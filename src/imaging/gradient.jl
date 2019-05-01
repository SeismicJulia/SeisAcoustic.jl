"""
   compute the source-side wave field as (p[it+1/2]-p[it-1/2]-f[it]) / k
"""
function apply_image_condition!(g::Vector{Tv}, p_adj::Vector{Tv}, wfd1::Wavefield,
         wfd2::Wavefield, params::ModelParams) where {Tv <: AbstractFloat}

    # number of elements
    N = params.nz * params.nx

    # constant
    # double = params.data_format(2.0)

    # compute gradient
    for i = 1 : N
        # dpdt = (wfd2.p[i] - wfd1.p[i]) / (params.rho[i] * params.vel[i] * params.vel[i])
        dpdt = (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        g[i] = g[i] + dpdt * p_adj[i]
    end
    return nothing
end

"""
   compute the gradient of bulk modulus
"""
function bulk_gradient(rec::Recordings, path_bnd::Ts, path_wfd::Ts,
         src::Source, params::ModelParams) where {Ts <: String}

    # model length
    N = params.nz * params.nx

    # allocate memory for computing adjoint wavefield
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp  = zeros(params.data_format, params.Nz * params.Nx)
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

    # ============================================
    # pre = zeros(params.data_format, N * params.nt)
    # idx_o = N*params.nt - N + 1
    # copyto!(pre, idx_o, wfd2.p, 1, N)
    # ============================================

    # allocate memory to save the adjoint pressure field
    p_adj = zeros(params.data_format, N)

    # allocate memory to save the gradient
    g     = zeros(params.data_format, N)

    # inject the recordings to the last snapshot
    inject_rec2spt!(spt2, rec, params.nt)

    # get the adjoint wavefield at the last time step
    sample_adjoint2pre!(p_adj, spt2, params)

    # open the file of boundary wavefield
    fid_bnd = open(path_bnd, "r")
    bnd     = WavefieldBound(params)

    # prepare for backward reconstruction
    subtract_source!(wfd2, src, params.nt)
    read_one_boundary!(bnd, fid_bnd, params.nt-1, params)
    one_step_backward!(wfd1, wfd2, bnd, params,
                       wfd_z1, wfd_z2, wfd_x1, wfd_x2)

    # =================================
    # idx_o = idx_o - N
    # copyto!(pre, idx_o, wfd1.p, 1, N)
    # =================================

    # compute the gradient
    apply_image_condition!(g, p_adj, wfd1, wfd2, params)

    # prepare for the next step backward reconstruction
    copy_wavefield!(wfd2, wfd1)

    # backward propagation
    for it = params.nt-1 : -1 : 2

        # compute adjoint wavefield
        one_step_adjoint!(spt1, spt2, params,
                          tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, rec, it)
        copy_snapshot!(spt2, spt1)

        # get the adjoint pressure at the current time step
        sample_adjoint2pre!(p_adj, spt2, params)

        # read the source-side wavefield at the last time step
        subtract_source!(wfd2, src, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)
        one_step_backward!(wfd1, wfd2, bnd, params,
                           wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # =================================
        # idx_o = idx_o - N
        # copyto!(pre, idx_o, wfd1.p, 1, N)
        # =================================

        # compute the gradient
        apply_image_condition!(g, p_adj, wfd1, wfd2, params)

        # prepare for the next step backward reconstruction
        copy_wavefield!(wfd2, wfd1)
    end

    close(fid_bnd)
    return g
end


# prove
