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

    # just a single source
    if typeof(src) <: Source && eltype(irz) <: Real

       # get the path for get_reflections
       path_obs = joinpath(dir_obs, "reflection_1.bin")

       # wrap argument
       argument_collection = (path_obs=path_obs, receiver_z=irz, receiver_x=irx, location_flag=location_flag,
                              source=src, fidiff_hete=params_hete, fidiff_homo=params_homo)

       wrap_get_reflections(argument_collection)

    # multiple source parallel
    elseif typeof(src) == Vector{Source}

       # number of simulation
       ns = length(src)

       # wrap parameters into vector of named tuple
       argument_collection = Vector{NamedTuple}(undef, ns)
       for i = 1 : ns

           file_name = join(["reflection_" "$i" ".bin"])
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

        return nothing
    end

    # just a single source
    if typeof(src) <: Source

       # set the path
       path_bnd = joinpath(dir_sourceside, "boundary_1.rsf" )
       path_wfd = joinpath(dir_sourceside, "wavefield_1.rsf")
       path_sws = joinpath(dir_sourceside, "strength_1.rsf" )

       # wrap argument
       argument_collection = (source=src, fidiff=fidiff, path_bnd=path_bnd, path_wfd=path_wfd, path_sws=path_sws)
       wrap_get_wavefield_bound(argument_collection)

       # compute normalization factor
       nal  = zeros(fidiff.data_format, fidiff.nz, fidiff.nx)
       (hdr, sws) = read_RSdata(path_sws)
       nal .= 1.0 ./ sws

       # remove the source strength file
       remove_flag && rm(path_sws)

    # multiple source parallel
    elseif typeof(src) == Vector{Source}

       # number of simulation
       ns = length(src)

       # wrap parameters into vector of named tuple
       argument_collection = Vector{NamedTuple}(undef, ns)
       for i = 1 : ns

           file_name = join(["boundary_" "$i" ".rsf"])
           path_bnd  = joinpath(dir_sourceside, file_name)

           file_name = join(["wavefield_" "$i" ".rsf"])
           path_wfd  = joinpath(dir_sourceside, file_name)

           file_name = join(["strength_" "$i" ".rsf"])
           path_sws  = joinpath(dir_sourceside, file_name)

           # pack parameter for source-side wavefield
           argument_collection[i] = (source=src[i], fidiff=fidiff, path_bnd=path_bnd,
                                     path_wfd=path_wfd, path_sws=path_sws)
       end

       # do simulation parallel
       pmap(wrap_get_wavefield_bound, argument_collection)

       # compute normalization factor
       nal = zeros(fidiff.data_format, fidiff.nz, fidiff.nx)
       for i = 1 : ns
           (hdr, sws) = read_RSdata(argument_collection[i].path_sws)
           nal .= nal .+ sws

           # remove the source strength file
           remove_flag && rm(argument_collection[i].path_sws)
       end
       nal .= 1.0 ./ nal
    end

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
function get_born_forward(dir_born::Ts, irz::Ti, irx::Ti, m::Vector{Tv}, dir_bnd::Ts,
                          src::T, fidiff::TdParams; location_flag="index", path_precondition="NULL") where {Ts<:String, Ti<:Vector, Tv<:AbstractFloat, T<:Union{Source, Vector{Source}}}

    function wrap_born_forward(params::NamedTuple)

        rec = Recordings(params.receiver_z, params.receiver_x, params.fidiff; location_flag=params.location_flag)
        born_approximation_forward!(rec, params.m, params.path_bnd, params.source, params.fidiff)
        write_recordings(params.path_born, rec)

        return nothing
    end

    # apply preconditioner to model parameter
    p = copy(m)
    if path_precondition != "NULL"
       (hdr, scale) = read_RSdata(path_precondition)
        scale = vec(scale)
        p .= scale .* m
    end

    # just a single source
    if typeof(src) <: Source && eltype(irz) <: Real

       path_born = joinpath(dir_born, "born_1.bin")
       path_bnd  = joinpath(dir_bnd , "boundary_1.rsf")
       argument_collecttion = (path_born=path_born, receiver_z=irz, receiver_x=irx, location_flag=location_flag, m=p,
                               path_bnd=path_bnd, source=src, fidiff=fidiff)
       wrap_born_forward(argument_collection)

    # multiple sources parallel
    elseif typeof(src) == Vector{Source}

       ns = length(src)
       argument_collection = Vector{NamedTuple}(undef, ns)

       for i = 1 : ns

           file_name = join(["born_" "$i" ".bin"])
           path_born = joinpath(dir_born, file_name)

           file_name = join(["boundary_" "$i" ".rsf"])
           path_bnd  = joinpath(dir_bnd, file_name)

           # OBN acquisition geometry
           if eltype(irz) <: Real
              argument_collection[i] = (path_born=path_born, receiver_z=irz, receiver_x=irx, location_flag=location_flag, m=p,
                                        path_bnd=path_bnd, source=src[i], fidiff=fidiff)

           # towed streamer
           elseif eltype(irz) <: Vector
              argument_collection[i] = (path_born=path_born, receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag, m=p,
                                        path_bnd=path_bnd, source=src[i], fidiff=fidiff)
           else
              error("wrong type receiver locations")
           end
       end

       # do simulation parallel
       pmap(wrap_born_forward, argument_collection)

    end

    return nothing
end

"""
   the adjoint operator of born approximation
"""
function get_born_adjoint(path_m::Ts, dir_rec::Ts, dir_bnd::Ts, dir_wfd::Ts, src::T, fidiff::TdParams;
                          path_precondition="NULL", remove_flag=true) where {Ts<:String, T<:Union{Source, Vector{Source}}}

    function wrap_born_adjoint(params::NamedTuple)

        rec = read_recordings(params.path_rec)
        m   = born_approximation_adjoint(rec, params.path_bnd, params.path_wfd, params.source, params.fidiff)
        hdr = RegularSampleHeader(m; title="temporary file")
        write_RSdata(params.path_tmp, hdr, m)

        return nothing
    end

    # the directory to save model
    dir_tmp = splitdir(path_m)[1]

    # just a single source
    if typeof(src) <: Source && eltype(irz) <: Real

       path_tmp = joinpath(dir_tmp, "temporary_1.rsf")
       path_rec = joinpath(dir_rec, "recording_1.bin")
       path_bnd = joinpath(dir_bnd, "boundary_1.rsf" )
       path_wfd = joinpath(dir_wfd, "wavefield_1.rsf")

       argument_collecttion = (path_tmp=path_tmp, path_rec=path_rec, path_bnd=path_bnd,
                               path_wfd=path_wfd, source=src, fidiff=fidiff)
       wrap_born_adjoint(argument_collection)

       # read the adjoint result
       (hdr, m) = read_RSdata(path_tmp)

       # remove the temporary file
       remove_flag && rm(argument_collection.path_tmp)

    # multiple sources parallel
    elseif typeof(src) == Vector{Source}

       ns = length(src)
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

           argument_collecttion[i] = (path_tmp=path_tmp, path_rec=path_rec, path_bnd=path_bnd,
                                      path_wfd=path_wfd, source=src[i], fidiff=fidiff)
       end

       # do simulation parallel
       pmap(wrap_born_adjoint, argument_collection)

       # read the adjoint result
       (hdr, p) = read_RSdata(argument_collection[1].path_tmp)
       remove_flag && rm(argument_collection[1].path_tmp)

       for i = 2 : ns
           (hdr, tmp) = read_RSdata(argument_collection[i].path_tmp)
           rm(argument_collection[i].path_tmp)
           p .= p .+ tmp
       end
    end

    # apply preconditioner to model parameter
    if path_precondition != "NULL"
       (hdr, scale) = read_RSdata(path_precondition)
       scale = vec(scale)
       p   .*= scale
    end

    hdr = RegularSampleHeader(p, title="image")
    write_RSdata(path_m, hdr, m)

    return nothing

end
