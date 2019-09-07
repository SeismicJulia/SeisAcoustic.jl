"""
   get the reflections for multiple single source, parallel computation is implemented
"""
function get_reflections(path::Ts, irz::Ti, irx::Ti, src::T, params_hete::P, params_homo::P,
                         location_flag="index") where {Ts<:String, Ti<:Vector, T<:Union{Source, Vector{Source}}, P<:TdParams}

    length(irz) == length(irx) || error("number of receiver doesn't match")

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
        write_recordings(params.path_out, dobs)
    end

    # just a single source
    if typeof(src) <: Source && eltype(irz) <: Real

       # get the path for get_reflections
       path_refl = joinpath(path, "reflections_1.bin")

       # wrap argument
       argument_collection = (receiver_z=irz, receiver_x=irx, location_flag=location_flag, source=src,
                              fidiff_hete=params_hete, fidiff_homo=params_homo, path_out=path_refl)

       wrap_get_reflections(argument_collection)

    # multiple source parallel
  elseif typeof(src) == Vector{Source}

       # number of simulation
       ns = length(src)

       # wrap parameters into vector of named tuple
       argument_collection = Vector{NamedTuple}(undef, ns)
       for i = 1 : ns

           file_name = join(["reflections_" "$i" ".bin"])
           path_refl = joinpath(path, file_name)

           # OBN acquisition geometry
           if eltype(irz) <: Real
              argument_collection[i] = (receiver_z=irz, receiver_x=irx, location_flag=location_flag, source=src[i],
                                        fidiff_hete=params_hete, fidiff_homo=params_homo, path_out=path_refl)

           # towed streamer
           elseif eltype(irz) <: Vector
              argument_collection[i] = (receiver_z=irz[i], receiver_x=irx[i], location_flag=location_flag, source=src[i],
                                        fidiff_hete=params_hete, fidiff_homo=params_homo, path_out=path_refl)
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
function get_wavefield_bound(path::Ts, src::T, fidiff::P; remove_flag=true) where {Ts<:String, T<:Union{Source, Vector{Source}}, P<:TdParams}

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
       path_bnd = joinpath(path, "boundary_1.rsf"   )
       path_wfd = joinpath(path, "wavefield_1.rsf"  )
       path_sws = joinpath(path, "strength_1.rsf"   )
       path_nal = joinpath(path, "normalization.rsf")

       # wrap argument
       argument_collection = (source=src, fidiff=fidiff, path_bnd=path_bnd, path_wfd=path_wfd, path_sws=path_sws)
       wrap_get_wavefield_bound(argument_collection)

       # compute the normalization factor
       (hdr, sws) = read_RSdata(path_sws)

       nal  = zeros(fidiff.data_format, fidiff.nz, fidiff.nx)
       nal .= 1.0 ./ sws
       hdr  = RegularSampleHeader(nal; title="source normalization")
       write_RSdata(path_nal, hdr, nal)

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
           path_bnd  = joinpath(path, file_name)

           file_name = join(["wavefield_" "$i" ".rsf"])
           path_wfd  = joinpath(path, file_name)

           file_name = join(["strength_" "$i" ".rsf"])
           path_sws  = joinpath(path, file_name)

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

       # save
       hdr = RegularSampleHeader(nal; title="source normalization")
       path_nal = joinpath(path, "normalization.rsf")
       write_RSdata(path_nal, hdr, nal)
    end

    return nothing
end
