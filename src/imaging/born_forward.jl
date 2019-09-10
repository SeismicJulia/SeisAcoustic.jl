"""
   Add the virtual source to pressure (just add to pz), the source-side wavefield
of all time-steps are saved in memory.
"""
function add_virtual_source!(spt::Snapshot, m::Vector{Tv}, wfd1::Wavefield, wfd2::Wavefield,
                             params::TdParams) where {Tv<:AbstractFloat}

    # total number of element
    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]
        dpdt = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        spt.px[j] = spt.px[j] +  dpdt * m[i]
    end

    return nothing
end

"""
   reconstruct source-side wavefield forwardly, the output can be used for born
approximation
"""
function sourceside_reconstruct_forward(path_bnd::Ts, src::Source,
         params::TdParams) where {Ts <: String}

    # length of one-step pressure field
    N = params.nz * params.nx

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = Wavefield(params)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # save the sourceside wavefield at one time step
    dpdt   = zeros(params.data_format, N)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd = WavefieldBound(params)

    # allocate memory for saving pressure field
    pre = zeros(params.data_format, N * (params.nt-1))

    # add source to the first wavefield
    add_source!(wfd1, src, 1)

    # the start index for saving the current source-side wavefield
    idx_o = 1

    # loop over time stepping
    for it = 2 : params.nt

        # read the boundary value
        read_one_boundary!(bnd, fid_bnd, it, params)

        # forward time steping and correcting boundaries
        one_step_forward!(wfd2, wfd1, bnd, params,
                          tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # subtract the source term if source located in the saved boundary area
        if (src.isz <= params.order || src.isz > params.nz-params.order ||
            src.isx <= params.order || src.isx > params.nx-params.order)

           subtract_source!(wfd2, src, it)
        end

        # compute the source-side wave field
        for i = 1 : N
            dpdt[i] = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        end

        # save the source-side field
        copyto!(pre, idx_o, dpdt, 1, N)
        idx_o = idx_o + N

        # add source to wavefield and prepare for next time step
        add_source!(wfd2, src, it)

        # prepare for next step
        copy_wavefield!(wfd1, wfd2)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt-1)

end

function sourceside_reconstruct_forward(path_bnd::Ts, srcs::Vector{Source},
         params::TdParams) where {Ts <: String}

    # length of one-step pressure field
    N = params.nz * params.nx

    # number of sources
    ns = length(srcs)

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = Wavefield(params)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # save the sourceside wavefield at one time step
    dpdt   = zeros(params.data_format, N)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd = WavefieldBound(params)

    # allocate memory for saving pressure field
    pre = zeros(params.data_format, N * (params.nt-1))

    # add source to the first wavefield
    add_multi_sources!(wfd1, srcs, 1)

    # the start index for saving the current source-side wavefield
    idx_o = 1

    # loop over time stepping
    for it = 2 : params.nt

        # read the boundary value
        read_one_boundary!(bnd, fid_bnd, it, params)

        # forward time steping and correcting boundaries
        one_step_forward!(wfd2, wfd1, bnd, params,
                          tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # subtract the source term if source located in the saved boundary area
        for i = 1 : ns

            if (srcs[i].isz <= params.order || srcs[i].isz > params.nz-params.order ||
                srcs[i].isx <= params.order || srcs[i].isx > params.nx-params.order)

                subtract_source!(wfd2, srcs[i], it)
            end
        end

        # compute the source-side wave field
        for i = 1 : N
            dpdt[i] = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        end

        # save the source-side field
        copyto!(pre, idx_o, dpdt, 1, N)
        idx_o = idx_o + N

        # add source to wavefield and prepare for next time step
        add_multi_sources!(wfd2, srcs, it)

        # prepare for next step
        copy_wavefield!(wfd1, wfd2)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt-1)

end

"""
   the first sample of recordings can be anything, the scattered wavefield is computed
from the second time step.
"""
function born_approximation_forward!(rec::Recordings, m::Vector{Tv}, path_bnd::Ts,
                                     src::Source, params::TdParams) where {Ts<:String, Tv<:AbstractFloat}

    # model length
    N = params.nz * params.nx

    # allocate memory for computing adjoint wavefield
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # allocate memory for reconstructing source-side wavefield backward
    wfd1 = Wavefield(params); add_source!(wfd1, src, 1); # get the wavefield at first time step
    wfd2 = Wavefield(params)
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # forward time stepping
    for it = 2 : params.nt

        # time stepping of scattered wavefield
        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # reconstruct source-side wavefield
        read_one_boundary!(bnd, fid_bnd, it, params)
        one_step_forward!(wfd2, wfd1, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # subtract source term if it being added during update boundary part
        if (src.isz <= params.order || src.isz > params.nz-params.order ||
            src.isx <= params.order || src.isx > params.nx-params.order)

            subtract_source!(wfd2, src, it)
        end

        # add the virtual source
        add_virtual_source!(spt2, m, wfd1, wfd2, params)

        # sampling scatter wavefield
        sample_spt2rec!(rec, spt2, it)

        # prepare for the next step backward reconstruction
        add_source!(wfd2, src, it)
        copy_wavefield!(wfd1, wfd2)
        copy_snapshot!(spt1, spt2)

    end

    # close the boundary value file
    close(fid_bnd)

    return nothing
end

"""
   Born approximation for sourceside wavefield is generated by simultaneouse source
"""
function born_approximation_forward!(rec::Recordings, m::Vector{Tv}, path_bnd::Ts,
                                    srcs::Vector{Source}, params::TdParams) where {Ts<:String, Tv<:AbstractFloat}

    # model length
    N = params.nz * params.nx

    # number of sources
    ns= length(srcs)

    # allocate memory for computing adjoint wavefield
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # allocate memory for reconstructing source-side wavefield backward
    wfd1 = Wavefield(params); add_multi_sources!(wfd1, srcs, 1); # get the wavefield at first time step
    wfd2 = Wavefield(params)
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # forward time stepping
    for it = 2 : params.nt

        # time stepping of scattered wavefield
        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # reconstruct source-side wavefield
        read_one_boundary!(bnd, fid_bnd, it, params)
        one_step_forward!(wfd2, wfd1, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # subtract source term if it being added during update boundary part
        # loop over all the sources
        for i = 1 : ns

            if (srcs[i].isz <= params.order || srcs[i].isz > params.nz-params.order ||
                srcs[i].isx <= params.order || srcs[i].isx > params.nx-params.order)

                subtract_source!(wfd2, srcs[i], it)
            end
        end


        # add the virtual source
        add_virtual_source!(spt2, m, wfd1, wfd2, params)

        # sampling scatter wavefield
        sample_spt2rec!(rec, spt2, it)

        # prepare for the next step backward reconstruction
        add_multi_sources!(wfd2, srcs, it)
        copy_wavefield!(wfd1, wfd2)
        copy_snapshot!(spt1, spt2)

    end

    # close the boundary value file
    close(fid_bnd)

    return nothing
end
