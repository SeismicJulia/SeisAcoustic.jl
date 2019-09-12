function born_approximation(output::Ts, input::Ts, iflag::Ti;
         irz::Vector, irx::Vector, location_flag="index", dir_sourceside="NULL", normalization_flag=true, mute_index=0, fidiff::TdParams) where {Ts<:String, Ti<:Int64}

    # clean the space
    rm(output, force=true, recursive=true)

    if iflag == 1
       born_approximation_forward!(output, input, irz, irx, dir_sourceside, fidiff;
                                  location_flag=location_flag, normalization_flag=normalization_flag, mute_index=mute_index)

    elseif iflag == 2
       born_approximation_adjoint(output, input, dir_sourceside, fidiff;
                                  normalization_flag=normalization_flag, mute_index=mute_index)

    else
       error("undefined operator behavior")
    end

    return nothing
end

"""
   perform y = ax + by for multiple recordings
"""
function recordings_axpby!(a, dir_x::Ts, b, dir_y::Ts) where {Ts<:String}

    # determine number of files
    file_name = readdir(dir_x)
    nx        = length(file_name)
    if file_name[1] == ".DS_Store" # in case of MacOS
       nx = nx - 1
    end

    file_name = readdir(dir_y)
    ny        = length(file_name)
    if file_name[1] == ".DS_Store" # in case of MacOS
       ny = ny - 1
    end

    nx == ny || error("number of files doesn't match")

    r = 0.0
    for i = 1 : nx
        file_name = join(["recordings_" "$i" ".bin"])

        path_x    = joinpath(dir_x, file_name)
        path_y    = joinpath(dir_y, file_name)

        rec_x     = read_recordings(path_x)
        rec_y     = read_recordings(path_y)

        for i2 = 1 : rec_x.nr
            for i1 = 1 : rec_x.nt
                rec_y.p[i1,i2] = a*rec_x.p[i1,i2] + b*rec_y.p[i1,i2]
            end
        end

        write_recordings(path_y, rec_y)
    end

    return nothing
end

"""
   perform y = ax + by for variables in model space
"""
function image_axpby!(a, path_x::Ts, b, path_y::Ts) where {Ts<:String}

    (hdr1, x) = read_RSdata(path_x)
    (hdr2, y) = read_RSdata(path_y)

    for i = 1 : length(x)
        y[i] = a*x[i] + b*y[i]
    end

    write_RSdata(path_y, hdr2, y)

    return nothing
end

function recordings_norm(dir::Ts) where {Tv<:Number, Ts<:String}

    # determine number of files
    file_name = readdir(dir)
    n         = length(file_name)
    if file_name[1] == ".DS_Store" # in case of MacOS
       n = n - 1
    end

    r = 0.0
    for i = 1 : n
        file_name = join(["recordings_" "$i" ".bin"])
        path      = joinpath(dir, file_name)
        rec       = read_recordings(path)
        tmp       = norm(rec.p)
        r        += tmp * tmp
    end

    return sqrt(r)
end
