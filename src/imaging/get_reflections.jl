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

           file_name = join(["reflections_" "$i" ".rec"])
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

end
