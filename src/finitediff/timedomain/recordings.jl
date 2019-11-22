"""
   structure for recordings, only pressure field are sampled
"""
struct Recordings{Ti<:Int64, Tv<:AbstractFloat}
     nt      :: Ti          # number of samples per trace
     nr      :: Ti          # number of receiver
     dt      :: Tv
     irz     :: Vector{Ti}  # vertical grid index of receivers
     irx     :: Vector{Ti}  # horizontal grid index of receivers
     spt2rec :: Vector{Ti}  # index mapping of receiver to snapshot
     p       :: Matrix{Tv}  # recordings of pressure field
end

"""
   constructor for recordings
"""
function Recordings(rz::Vector, rx::Vector,
         params::TdParams; location_flag="index")

    # number of receivers
    nr  = length(rz)
    nr == length(rx) || error("length of receiver coordinate doesn't match")

    irz = zeros(Int64, nr)
    irx = zeros(Int64, nr)

    # location are provided as index
    if location_flag == "index"
       for i = 1 : nr
           irz[i] = round(Int64, rz[i])
           irx[i] = round(Int64, rx[i])
       end

    # location are given as distance
    elseif location_flag == "distance"
       for i = 1 : nr
           irz[i] = round(Int64, rz[i]/params.dz) + 1
           irx[i] = round(Int64, rx[i]/params.dx) + 1
       end

    else
       error("wrong specification of receiver location")
    end

    # error checking
    for i = 1 : nr
        if irz[i] > params.nz || irx[i] > params.nx || irz[i] < 1 || irx[i] < 1
           error("receiver located outside of modeling area")
        end
    end

    # can't put receivers on the free surface
    if params.free_surface
       for i = 1 : nr
           if irz[i] == 1
              error("can't put receiver at free surface, move it deeper")
           end
       end
    end

    # the auxillary vector mapping snapshot to recordings
    spt2rec = zeros(Int64, nr)
    for i = 1 : nr
        spt2rec[i] = (irx[i]+params.npml-1) * params.Nz + irz[i] + params.ntop
    end

    # initalize empty recordings
    p  = zeros(params.data_format, params.nt, nr)
    rec= Recordings(params.nt, nr, params.dt, irz, irx, spt2rec, p)

    return rec
end

"""
   write recordings into su file
"""
function write_recordings(path::String, rec::Recordings)

    fid = open(path, "w")

    # write number of samples
    write(fid, rec.nt)

    # write number of receivers
    write(fid, rec.nr)

    # write the locations of receivers (index)
    write(fid, rec.irz)
    write(fid, rec.irx)

    # write the index mapping vector
    write(fid, rec.spt2rec)

    # data format code
    if eltype(rec.p) == Float32
       write(fid, 1)             # 1 indicate Float32

    elseif eltype(rec.p) == Float64
       write(fid, 2)             # 2 indicate Float64

    else
       error("non-support data format")
    end

    # write time sampling interval
    write(fid, rec.dt)

    # write the recordings
    write(fid, rec.p)

    #close file
    close(fid)

    return nothing
end

"""
   read one shot gather
"""
function read_recordings(path::String)

    # open the file
    fid = open(path, "r")

    # number of samples
    nt  = read(fid, Int64)

    # number of traces
    nr  = read(fid, Int64)

    # read the index locations of receivers
    irz = zeros(Int64, nr); read!(fid, irz)
    irx = zeros(Int64, nr); read!(fid, irx)

    # read the index mapping vector
    spt2rec = zeros(Int64, nr); read!(fid, spt2rec)

    # read the data_format code
    code = read(fid, Int64)
    if code == 1
       data_format = Float32
    elseif code == 2
       data_format = Float64
    else
       error("non-support data format")
    end

    # if 8*(3*nr+2) + 8*(nr*nr+1) > filesize(path)
    #    data_format = Float32
    # else
    #    data_format = Float64
    # end

    # read the time sampling interval
    dt = read(fid, data_format)

    # read the recordings
    p = zeros(data_format, nt * nr)
    read!(fid, p); p = reshape(p, nt, nr);

    #close file
    close(fid)

    # organize them into a structure
    rec = Recordings(nt, nr, dt, irz, irx, spt2rec, p)

    return rec
end

"""
   sampling one snashot to fill the recordings
"""
function sample_spt2rec!(rec::Recordings, spt::Snapshot, it::Int64)

    # loop over receivers
    for i = 1 : rec.nr
        idx = rec.spt2rec[i]
        rec.p[it,i] = spt.pz[idx] + spt.px[idx]
    end

    return nothing
end

"""
   the adjoint of sampling operator, inject recordings to snapshot
"""
function inject_rec2spt!(spt::Snapshot, rec::Recordings, it::Int64)

    # loop over receivers
    for i = 1 : rec.nr
        idx = rec.spt2rec[i]
        spt.pz[idx] = spt.pz[idx] + rec.p[it,i]
        spt.px[idx] = spt.px[idx] + rec.p[it,i]
    end

    return nothing
end

"""
   check whether two shot gathers are equal to each other
"""
function recordings_isequal(rec1::Tr, rec2::Tr) where {Tr<:Recordings}

    rec1.nt      != rec2.nt      && (return false)
    rec1.nr      != rec2.nr      && (return false)
    rec1.dt      != rec2.dt      && (return false)
    rec1.irz     != rec2.irz     && (return false)
    rec1.irx     != rec2.irx     && (return false)
    rec1.spt2rec != rec2.spt2rec && (return false)
    norm(rec1.p-rec2.p)/norm(rec1.p) > 10*eps(eltype(rec1.p)) && (return false)

    return true
end
