"""
   compute the source-side wave field as (p[it+1/2]-p[it-1/2]-f[it]) / k
"""
function apply_image_condition!(g::Vector{Tv}, spt::Snapshot, wfd1::Wavefield,
         wfd2::Wavefield, params::TdParams) where {Tv <: AbstractFloat}

    # number of elements
    N = params.nz * params.nx

    # update gradient by current adjoint snapshot
    for i = 1 : N
        dpdt = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        j    = params.spt2wfd[i]
        g[i] = g[i] + dpdt * spt.px[j]
    end

    return nothing
end

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
   reconstruct source-side backward, which can be used for the adjoint of born
approximation
"""
function sourceside_reconstruct_backward(path_bnd::Ts, path_wfd::Ts,
         srcs::Vector{Source}, params::TdParams) where {Ts <: String}

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
        subtract_multi_sources!(wfd2, srcs, it+1)

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

"""
   the adjoint operator of born approximation
"""
function born_approximation_adjoint(path_m::Ts, dir_rec::Ts, dir_sourceside::Ts, fidiff::TdParams;
                          normalization_flag=true, remove_flag=true, mute_index::Int64=0) where {Ts<:String, T<:Union{Source, Vector{Source}}}

    function wrap_born_adjoint(params::NamedTuple)

        rec = read_recordings(params.path_rec)
        src = read_source(params.path_src)

        m   = born_approximation_adjoint(rec, params.path_bnd, params.path_wfd, src, params.fidiff)
        hdr = RegularSampleHeader(m; title="temporary file")
        write_RSdata(params.path_tmp, hdr, m)

        return nothing
    end

    # set up directory
    dir_bnd = joinpath(dir_sourceside, "boundary")
    dir_wfd = joinpath(dir_sourceside, "wavefield")
    dir_src = joinpath(dir_sourceside, "source")
    isdir(dir_rec) || error("reflections doesn't exist")
    isdir(dir_bnd) || error("boundary file doesn't exist")
    isdir(dir_wfd) || error("wavefield file doesn't exist")
    isdir(dir_src) || error("source wavelet file doesn't exist")

    # the directory to save model
    dir_m = splitdir(path_m)[1]
    if !isdir(dir_m)
       mkdir(dir_m)
       if !isdir(dir_m) # check the directory is created
          error("can't create directory for the variable in model space")
       end
    end

    dir_tmp = joinpath(dir_m, "temporary")
    if !isdir(dir_tmp)
       mkdir(dir_tmp)
       if !isdir(dir_tmp) # check the directory is created
          error("can't create directory for one shot adjoint born approximation")
       end
    end

    # determine the number of shot
    file_rec = readdir(dir_rec)
    ns       = length(file_rec)
    if file_rec[1] == ".DS_Store"
       ns = ns - 1
    end

    file_src = readdir(dir_src)
    nb       = length(file_src)
    if file_src[1] == ".DS_Store"
       nb = nb - 1
    end
    ns == nb || error("number of shot doesn't match")

    argument_collection = Vector{NamedTuple}(undef, ns)
    for i = 1 : ns

        file_name = join(["temporary_" "$i" ".bin"])
        path_tmp  = joinpath(dir_tmp, file_name)

        file_name = join(["recordings_" "$i" ".bin"])
        path_rec  = joinpath(dir_rec, file_name)

        file_name = join(["boundary_" "$i" ".rsf"])
        path_bnd  = joinpath(dir_bnd, file_name)

        file_name = join(["wavefield_" "$i" ".rsf"])
        path_wfd  = joinpath(dir_wfd, file_name)

        file_name = join(["source_" "$i" ".bin"])
        path_src  = joinpath(dir_src, file_name)

        argument_collection[i] = (path_tmp=path_tmp, path_rec=path_rec, path_bnd=path_bnd,
                                  path_wfd=path_wfd, path_src=path_src, fidiff=fidiff)
    end

    # do simulation parallel
    if nprocs() == 1
       for i = 1 : ns
           wrap_born_adjoint(argument_collection[i])
       end
    else
       pmap(wrap_born_adjoint, argument_collection)
    end

    # read the adjoint result
    (hdr, p) = read_RSdata(argument_collection[1].path_tmp)
    for i = 2 : ns
        (hdr, tmp) = read_RSdata(argument_collection[i].path_tmp)
        p .= p .+ tmp
    end

    # remove the temporary folder
    remove_flag && rm(dir_tmp, force=true, recursive=true)

    # apply muting
    if mute_index > 0
       p  = reshape(p, fidiff.nz, fidiff.nx)
       p[1:mute_index,:] .= 0.0
       p  = vec(p)
    end

    # apply preconditioner to model parameter
    if normalization_flag
       path_normalization = joinpath(dir_sourceside, "normalization.rsf")
       (hdr, scale) = read_RSdata(path_normalization)
       p .*= vec(scale)
    end

    hdr = RegularSampleHeader(p, title="image")
    write_RSdata(path_m, hdr, p)

    return nothing
end
