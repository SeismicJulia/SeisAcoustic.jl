"""
   get the observations for multiple single source, parallel computation is implemented,
the size of velocity model is changed according to the acquisition geometry
"""
function get_shotgather(dir_obs::Ts, isz, isx, wlet, irz::Ti, irx::Ti,
                        vel::Matrix{Tv}, rho::Matrix{Tv},  dz, dx, dt, tmax;
                        location_flag="index",  data_format=Float32, order=2, free_surface=true, npad=100) where {Ts<:String, Ti<:Vector, Tv<:Real}

    # get the original model size
    (nz, nx) = size(vel)
    (nz, nx) == size(rho) || error("model size doesn't match")

    # define function accept named tuple
    function wrap_get_shotgather(params::NamedTuple)

        # determine model range in spatial direction
        # lower end
        ix_lower = minimum(vcat(params.source_x, params.receiver_x))
        if ix_lower > params.npad
           ix_lower = ix_lower - params.npad
        else
           ix_lower = 1
        end
        # upper end
        ix_upper = maximum(vcat(params.source_x, params.receiver_x))
        if ix_upper + params.npad > nx
           ix_upper = nx
        else
           ix_upper = ix_upper + params.npad
        end

        # grab part of velocity model (crop only in horizontal direction)
        vel = params.vel[:, ix_lower:ix_upper]
        rho = params.rho[:, ix_lower:ix_upper]
        fidiff = TdParams(rho, vel, params.free_surface, params.dz, params.dx, params.dt, params.tmax;
                          data_format=params.data_format, order=params.order)

        # adjust receiver's location
        receiver_x = copy(params.receiver_x)
        for i = 1 : length(receiver_x)
            receiver_x[i] = receiver_x[i] - ix_lower + 1
        end
        dobs = Recordings(params.receiver_z, receiver_x, fidiff; location_flag=params.location_flag)

        # adjust source's location
        source_x = params.source_x - ix_lower + 1
        src = Source(params.source_z, source_x, fidiff; p=params.wlet)

        # do simulation
        multi_step_forward!(dobs, src, fidiff)

        # save the result to disk
        write_recordings(params.path_obs, dobs)
    end

    # create a fold save observations
    if isdir(dir_obs)
       rm(dir_obs, force=true, recursive=true)
    end
    mkdir(dir_obs)
    if !isdir(dir_obs)
       error("can't create directory for reflections")
    end

    # determine number of shot
    ns =  length(isx)
    ns == length(isz) || error("source coordinate doesn't match")

    # wrap parameters into vector of named tuple
    argument_collection = Vector{NamedTuple}(undef, ns)
    for i = 1 : ns

        file_name = join(["recordings_" "$i" ".bin"])
        path_obs  = joinpath(dir_obs, file_name)

        # OBN acquisition geometry
        if eltype(irz) <: Real
           argument_collection[i] = (path_obs=path_obs, receiver_z=irz, receiver_x=irx, location_flag=location_flag, source_z=isz[i], source_x=isx[i], wlet=wlet,
                                     free_surface=free_surface, dz=dz, dx=dx, dt=dt, tmax=tmax, data_format=data_format, order=order, vel=vel, rho=rho, npad=npad)

        # towed streamer
        elseif eltype(irz) <: Vector
           argument_collection[i] = (path_obs=path_obs, receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag, source_z=isz[i], source_x=isx[i], wlet=wlet,
                                     free_surface=free_surface, dz=dz, dx=dx, dt=dt, tmax=tmax, data_format=data_format, order=order, vel=vel, rho=rho, npad=npad)
        else
           error("wrong type receiver locations")
        end
    end

    # do simulation parallel
    if nprocs() == 1
       for i = 1 : ns
           wrap_get_shotgather(argument_collection[i])
       end
    else
       pmap(wrap_get_shotgather, argument_collection)
    end

    return nothing
end

"""
   get the observations for multiple single source, parallel computation is implemented
"""
function get_observations(dir_obs::Ts, irz::Ti, irx::Ti, src::T, fidiff::P,
                         location_flag="index") where {Ts<:String, Ti<:Vector, T<:Union{Source, Vector{Source}}, P<:TdParams}

    # define function accept named tuple
    function wrap_get_observations(params::NamedTuple)

        # get whole observations
        dobs = Recordings(params.receiver_z, params.receiver_x, params.fidiff; location_flag=params.location_flag)

        # do simulation
        multi_step_forward!(dobs, params.source, params.fidiff)

        # save the result to disk
        write_recordings(params.path_obs, dobs)
    end

    # create a fold save observations
    if isdir(dir_obs)
       rm(dir_obs, force=true, recursive=true)
    end
    mkdir(dir_obs)
    if !isdir(dir_obs)
       error("can't create directory for reflections")
    end

    # determine number of shot
    if typeof(src) <: Source
       ns  = 1
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
                                     source=src[i], fidiff=fidiff)

        # towed streamer
        elseif eltype(irz) <: Vector
           argument_collection[i] = (path_obs=path_obs, receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag,
                                     source=src[i], fidiff=fidiff)
        else
           error("wrong type receiver locations")
        end
    end

    # do simulation parallel
    if nprocs() == 1
       for i = 1 : ns
           wrap_get_observations(argument_collection[i])
       end
    else
       pmap(wrap_get_observations, argument_collection)
    end

    return nothing
end

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

    # create a fold save observations
    rm(dir_obs, force=true, recursive=true)
    mkdir(dir_obs)
    if !isdir(dir_obs)
       error("can't create directory for reflections")
    end

    # determine number of shot
    if typeof(src) <: Source
       ns  = 1
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
           argument_collection[i] = (path_obs=path_obs, receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag,
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

    # create a folder for sourceside wavefield
    rm(dir_sourceside, force=true, recursive=true)
    mkdir(dir_sourceside)
    if !isdir(dir_sourceside) # check the directory is created
       error("can't create directory for source-side wavefield")
    end

    # create folders to save the file
    dir_bnd = joinpath(dir_sourceside, "boundary")
    dir_wfd = joinpath(dir_sourceside, "wavefield")
    dir_sws = joinpath(dir_sourceside, "strength")
    dir_src = joinpath(dir_sourceside, "source")

    # if this folder doesn't exist, create one
    mkdir(dir_bnd)
    if !isdir(dir_bnd) # check the directory is created
       error("can't create directory for boundary")
    end

    mkdir(dir_wfd)
    if !isdir(dir_wfd) # check the directory is created
       error("can't create directory for wavefield")
    end

    mkdir(dir_sws)
    if !isdir(dir_sws) # check the directory is created
       error("can't create directory for source strength")
    end

    mkdir(dir_src)
    if !isdir(dir_src) # check the directory is created
       error("can't create directory for source wavelet")
    end

    # determine number of shot
    if typeof(src) <: Source
       ns  = 1
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
   create a project folder for LSRTM
"""
function initialize_lsrtm(work_dir::String)

    # check physical model
    dir_physical = joinpath(work_dir, "physical_model")
    isdir(dir_physical) || error("no physical model")

    dir_born = joinpath(work_dir, "data_space")
    rm(dir_born, force=true, recursive=true)
    mkdir(dir_born)
    if !isdir(dir_born) # check the directory is created
       error("can't create directory for variables in data space")
    end

    dir_model = joinpath(work_dir, "model_space")
    rm(dir_model, force=true, recursive=true)
    mkdir(dir_model)
    if !isdir(dir_model) # check the directory is created
       error("can't create directory for variables in model space")
    end

    # dir_sourceside = joinpath(work_dir, "sourceside")
    # rm(dir_sourceside, force=true, recursive=true)
    # mkdir(dir_sourceside)
    # if !isdir(dir_sourceside) # check the directory is created
    #    error("can't create directory for boundary values")
    # end

    return nothing
end
