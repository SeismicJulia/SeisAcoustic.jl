"""
   get the reflections for multiple single source, parallel computation is implemented
"""
function get_reflections(dir_obs::Ts, irz::Ti, irx::Ti, src::T, params_hete::P, params_homo::P;
                         location_flag="index") where {Ts<:String, Ti<:Vector, T<:Union{Source, Vector{Source}}, P<:TdParams}

    # define function accept named tuple
    function wrap_get_reflections(params::NamedTuple)

        # get whole observations
        dobs = Recordings(params.receiver_z, params.receiver_x, params.fidiff_hete; location_flag=params.location_flag)
        dire = Recordings(params.receiver_z, params.receiver_x, params.fidiff_homo; location_flag=params.location_flag)

        # do simulation
        multi_step_forward!(dobs, params.source, params.fidiff_hete)
        multi_step_forward!(dire, params.source, params.fidiff_homo)

        dobs.p .= dobs.p .- dire.p

        # save the result to disk
        write_recordings(params.path_obs, dobs)
    end

    # create a fold if it doesn't exist
    if !isdir(dir_obs)
       mkdir(dir_obs)
       if !isdir(dir_obs) # check the directory is created
          error("can't create directory for reflections")
       end
    end

    # determine number of shot
    if typeof(src) <: Source
       ns = 1
       src = [src]
    else
       ns = length(src)
    end

    # wrap parameters into vector of named tuple
    argument_collection = Vector{NamedTuple}(undef, ns)
    for i = 1 : ns

        file_name = join(["recordings_" "$i" ".bin"])
        path_obs  = joinpath(dir_obs, file_name)

        # OBN acquisition geometry
        if eltype(irz) <: Real
           argument_collection[i] = (path_obs=path_obs, receiver_z=irz, receiver_x=irx, location_flag=location_flag,
                                     source=src[i], fidiff_hete=params_hete, fidiff_homo=params_homo)

        # towed streamer
        elseif eltype(irz) <: Vector
           argument_collection[i] = (path_refl=path_refl, receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag,
                                     source=src[i], fidiff_hete=params_hete, fidiff_homo=params_homo)
        else
           error("wrong type receiver locations")
        end
    end

    # do simulation parallel
    if nprocs() == 1
       for i = 1 : ns
           wrap_get_reflections(argument_collection[i])
       end
    else
       pmap(wrap_get_reflections, argument_collection)
    end

    return nothing
end

"""
   get wavefield boundary, last wavefield and sourceside wavefield strength
"""
function get_wavefield_bound(dir_sourceside::Ts, src::T, fidiff::P; remove_flag=true) where {Ts<:String, T<:Union{Source, Vector{Source}}, P<:TdParams}

    # define function accept named tuple
    function wrap_get_wavefield_bound(params::NamedTuple)

        # do simulation
        multi_step_forward!(params.source, params.fidiff;
                            path_bnd=params.path_bnd, path_wfd=params.path_wfd, path_sws=params.path_sws)

        # save the source wavelet
        write_source(params.path_src, params.source)
        return nothing
    end

    # if this folder doesn't exist, create one
    if !isdir(dir_sourceside)
       mkdir(dir_sourceside)
       if !isdir(dir_sourceside) # check the directory is created
          error("can't create directory for source-side wavefield")
       end
    end

    # create folders to save the file
    dir_bnd = joinpath(dir_sourceside, "boundary")
    dir_wfd = joinpath(dir_sourceside, "wavefield")
    dir_sws = joinpath(dir_sourceside, "strength")
    dir_src = joinpath(dir_sourceside, "source")

    # if this folder doesn't exist, create one
    if !isdir(dir_bnd)
       mkdir(dir_bnd)
       if !isdir(dir_bnd) # check the directory is created
          error("can't create directory for boundary")
       end
    end

    if !isdir(dir_wfd)
       mkdir(dir_wfd)
       if !isdir(dir_wfd) # check the directory is created
          error("can't create directory for wavefield")
       end
    end

    if !isdir(dir_sws)
       mkdir(dir_sws)
       if !isdir(dir_sws) # check the directory is created
          error("can't create directory for source strength")
       end
    end

    if !isdir(dir_src)
       mkdir(dir_src)
       if !isdir(dir_src) # check the directory is created
          error("can't create directory for source wavelet")
       end
    end

    # determine number of shot
    if typeof(src) <: Source
       ns = 1
       src = [src]
    else
       ns = length(src)
    end

    # wrap parameters into vector of named tuple
    argument_collection = Vector{NamedTuple}(undef, ns)
    for i = 1 : ns

        file_name = join(["boundary_" "$i" ".rsf"])
        path_bnd  = joinpath(dir_bnd, file_name)

        file_name = join(["wavefield_" "$i" ".rsf"])
        path_wfd  = joinpath(dir_wfd, file_name)

        file_name = join(["strength_" "$i" ".rsf"])
        path_sws  = joinpath(dir_sws, file_name)

        file_name = join(["source_" "$i" ".bin"])
        path_src  = joinpath(dir_src, file_name)

        # pack parameter for source-side wavefield
        argument_collection[i] = (source=src[i], fidiff=fidiff, path_bnd=path_bnd,
                                  path_wfd=path_wfd, path_sws=path_sws, path_src=path_src)
    end

    # do simulation parallel
    if nprocs() == 1
       for i = 1 : ns
           wrap_get_wavefield_bound(argument_collection[i])
       end
    else
       pmap(wrap_get_wavefield_bound, argument_collection)
    end

    # compute normalization factor
    nal = zeros(fidiff.data_format, fidiff.nz, fidiff.nx)
    for i = 1 : ns
        (hdr, sws) = read_RSdata(argument_collection[i].path_sws)
        nal .= nal .+ sws
    end
    nal .= 1.0 ./ nal

    # remove the source strength file
    remove_flag && rm(dir_sws, force=true, recursive=true)

    # save the normalization factor
    hdr = RegularSampleHeader(nal; title="source normalization")
    path_nal = joinpath(dir_sourceside, "normalization.rsf")
    write_RSdata(path_nal, hdr, nal)

    return nothing
end

"""
   parallel born approximation for multiple sources, dir_born and dir_bnd are directory
which contain the born approximation and boundary files
"""
function get_born_forward(dir_born::Ts, irz::Ti, irx::Ti, m::Vector{Tv}, dir_sourceside::Ts,
                          fidiff::TdParams; location_flag="index", normalization_flag=true) where {Ts<:String, Ti<:Vector, Tv<:AbstractFloat, T<:Union{Source, Vector{Source}}}

    function wrap_born_forward(params::NamedTuple)

        # read source wavelet
        src = read_source(params.path_src)
        rec = Recordings(params.receiver_z, params.receiver_x, params.fidiff; location_flag=params.location_flag)

        born_approximation_forward!(rec, params.m, params.path_bnd, src, params.fidiff)
        write_recordings(params.path_born, rec)

        return nothing
    end

    # create folder to save the result
    if !isdir(dir_born)
       mkdir(dir_born)
       if !isdir(dir_born) # check the directory is created
          error("can't create directory for forward born approximation")
       end
    end

    dir_bnd = joinpath(dir_sourceside, "boundary")
    dir_src = joinpath(dir_sourceside, "source")
    isdir(dir_bnd) || error("boundary file doesn't exist")
    isdir(dir_src) || error("source wavelet file doesn't exist")

    # determine the number of shot
    file_bnd = readdir(dir_bnd)
    file_src = readdir(dir_src)
    ns       = length(file_src)
    if file_src[1] == ".DS_Store"
       ns = ns - 1
    end

    # apply preconditioner to model parameter
    p = copy(m)
    if normalization_flag
       path_normalization = joinpath(dir_sourceside, "normalization.rsf")
       (hdr, scale) = read_RSdata(path_normalization)
       scale = vec(scale)
       p .= scale .* m
    end

    # prepare argument
    argument_collection = Vector{NamedTuple}(undef, ns)
    for i = 1 : ns

        file_name = join(["recordings_" "$i" ".bin"])
        path_born = joinpath(dir_born, file_name)

        file_name = join(["boundary_" "$i" ".rsf"])
        path_bnd  = joinpath(dir_bnd, file_name)

        file_name = join(["source_" "$i" ".bin"])
        path_src  = joinpath(dir_src, file_name)

        # OBN acquisition geometry
        if eltype(irz) <: Real
           argument_collection[i] = (path_born=path_born, receiver_z=irz, receiver_x=irx, location_flag=location_flag, m=p,
                                     path_bnd=path_bnd, path_src=path_src, fidiff=fidiff)

        # towed streamer
        elseif eltype(irz) <: Vector
           argument_collection[i] = (path_born=path_born, receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag, m=p,
                                     path_bnd=path_bnd, path_src=path_src, fidiff=fidiff)
        else
           error("wrong type receiver locations")
        end
    end

    # do simulation parallel
    if nprocs() == 1
       for i = 1 : ns
           wrap_born_forward(argument_collection[i])
       end
    else
       pmap(wrap_born_forward, argument_collection)
    end

    return nothing
end

"""
   the adjoint operator of born approximation
"""
function get_born_adjoint(path_m::Ts, dir_rec::Ts, dir_sourceside::Ts, fidiff::TdParams;
                          normalization_flag=true, remove_flag=true) where {Ts<:String, T<:Union{Source, Vector{Source}}}

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

    file_bnd = readdir(dir_bnd)
    file_wfd = readdir(dir_wfd)
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

        file_name = join(["recording_" "$i" ".bin"])
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
    remove_flag && rm(dir_tmp, force-true, recursive=true)

    # apply preconditioner to model parameter
    if normalization_flag
       path_normalization = joinpath(dir_sourceside, "normalization.rsf")
       (hdr, scale) = read_RSdata(path_normalization)
       scale = vec(scale)
       p   .*= scale
    end

    hdr = RegularSampleHeader(p, title="image")
    write_RSdata(path_m, hdr, p)

    return nothing
end

"""
   create a project folder for LSRTM
"""
function initialize_lsrtm(work_dir::String)

    # check physical model
    dir_physical = joinpath(work_dir, "physical_model")
    isdir(dir_physical) || error("no physical model")

    dir_data = joinpath(work_dir, "data_space")
    rm(dir_data, force=true, recursive=true)
    mkdir(dir_data)
    if !isdir(dir_data) # check the directory is created
       error("can't create directory for variables in data space")
    end

    dir_model = joinpath(work_dir, "model_space")
    rm(dir_model, force=true, recursive=true)
    mkdir(dir_model)
    if !isdir(dir_model) # check the directory is created
       error("can't create directory for variables in model space")
    end

    dir_sourceside = joinpath(work_dir, "sourceside")
    rm(dir_sourceside, force=true, recursive=true)
    mkdir(dir_sourceside)
    if !isdir(dir_sourceside) # check the directory is created
       error("can't create directory for boundary values")
    end

    return nothing
end
