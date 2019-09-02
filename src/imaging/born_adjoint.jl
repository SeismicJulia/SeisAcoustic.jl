"""
   reconstruct source-side backward, which can be used for the adjoint of born
approximation
"""
function sourceside_reconstruct_backward(path_bnd::Ts, path_wfd::Ts,
         src::Source, params::TdParams) where {Ts <: String}

    # length of one-step pressure field
    N = params.nz * params.nx

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = read_one_wavefield(path_wfd, 1)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # memory for saving the sourceside wavefield at one time step
    dpdt   = zeros(params.data_format, N)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd = WavefieldBound(params)

    # pre-allocate memory for saving the reconstructed pressure field
    pre = zeros(params.data_format, N * (params.nt-1))
    idx_o = N * (params.nt-2) + 1

    # backward iterations
    for it = params.nt-1 : -1 : 1

        # remove dt * src[it+1]
        subtract_source!(wfd2, src, it+1)

        # read the boundary value
        read_one_boundary!(bnd, fid_bnd, it, params)

        # one step backward reconstruction
        one_step_backward!(wfd1, wfd2, bnd, params,
                           tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # compute the source-side wavefield
        for i = 1 : N
            dpdt[i] = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        end

        # save it
        copyto!(pre, idx_o, dpdt, 1, N)
        idx_o = idx_o - N

        # prepare for next step
        copy_wavefield!(wfd2, wfd1)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt-1)
end

"""
   the adjoint operator of born approximation of a single source
"""
function born_approximation_adjoint(rec::Recordings, path_bnd::Ts, path_wfd::Ts,
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

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # initialize gradient to be zeros
    m = zeros(params.data_format, N)

    # backward time stepping
    for it = params.nt : -1 : 2

        # compute adjoint wavefield
        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, rec, it)

        # reconstruct source-side wavefield
        subtract_source!(wfd2, src, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)
        one_step_backward!(wfd1, wfd2, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # compute the gradient
        apply_image_condition!(m, spt1, wfd1, wfd2, params)

        # prepare for the next step backward reconstruction
        copy_snapshot!(spt2, spt1)
        copy_wavefield!(wfd2, wfd1)

    end

    # close the boundary value file
    close(fid_bnd)

    # return the gradient of velocity model
    return m
end

"""
   the adjoint operator of born approximation of simultaneouse source
"""
function born_approximation_adjoint(rec::Recordings, path_bnd::Ts, path_wfd::Ts,
                           srcs::Vector{Source}, params::TdParams) where {Ts <: String}

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
    wfd2 = read_one_wavefield(path_wfd, 1)        # read the last wavefield
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # initialize gradient to be zeros
    m = zeros(params.data_format, N)

    # backward time stepping
    for it = params.nt : -1 : 2

        # compute adjoint wavefield
        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, rec, it)

        # reconstruct source-side wavefield
        subtract_multi_sources!(wfd2, srcs, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)
        one_step_backward!(wfd1, wfd2, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # compute the gradient
        apply_image_condition!(m, spt1, wfd1, wfd2, params)

        # prepare for the next step backward reconstruction
        copy_snapshot!(spt2, spt1)
        copy_wavefield!(wfd2, wfd1)

    end

    # close the boundary value file
    close(fid_bnd)

    # return the gradient of velocity model
    return m
end
