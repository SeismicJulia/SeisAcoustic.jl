"""
   Struct for snapshot which include four fields: vz, vx, pz, px
"""
struct Snapshot{Tv<:AbstractFloat}

    vz :: Vector{Tv}
    vx :: Vector{Tv}
    pz :: Vector{Tv}
    px :: Vector{Tv}
end

"""
   constructor of snapshot
"""
function Snapshot(params::TdParams)

    data_format = params.data_format
    N           = params.Nz * params.Nx

    return Snapshot(zeros(data_format,N), zeros(data_format,N),
                    zeros(data_format,N), zeros(data_format,N))
end

"""
   Struct for wavefiled which include three fields(vz, vx, p), PML boundary are truncated
"""
struct Wavefield{Tv<:AbstractFloat}
     vz :: Vector{Tv}
     vx :: Vector{Tv}
     p  :: Vector{Tv}
end

"""
   constructor of wavefield
"""
function Wavefield(params::TdParams)

    data_format = params.data_format
    N           = params.nz * params.nx

    return Wavefield(zeros(data_format,N), zeros(data_format,N), zeros(data_format,N))
end

"""
   sample snapshot to get wavefield
"""
function sample_spt2wfd(spt::Snapshot, params::TdParams)

    # initial a empty wavefield
    wfd = Wavefield(params)
    N   = length(params.spt2wfd)

    for i = 1 : N
        j = params.spt2wfd[i]
        wfd.vz[i] = spt.vz[j]
        wfd.vx[i] = spt.vx[j]
        wfd.p[i]  = spt.pz[j] + spt.px[j]
    end

    return wfd
end

"""
   sample snapshot to get pressure field
"""
function sample_adjoint2pre!(p::Vector{Tv}, spt::Snapshot,
         params::TdParams) where {Tv <: AbstractFloat}

    i = 0
    for j in params.spt2wfd
        i = i + 1
        p[i] = spt.px[j] # pz = px = p/2 inside the modeling area
    end
    return nothing
end

"""
   Create header of regularly sampled data for snapshot
"""
function snapshot_header(params::TdParams, interval::Int64)

    # the size of data
    n1 = params.Nz
    n2 = params.Nx
    n3 = 4       # include vz, vx, pz, px
    n4 = floor(Int64, (params.nt-1)/interval) + 1

    # label
    label1="z-axis"
    label2="x-axis"
    label3="component"
    label4="time step"

    return RegularSampleHeader(n1=n1, n2=n2, n3=n3, n4=n4,
           o1=0.0      , o2=0.0      , o3=1.0, o4=0.0    ,
           d1=params.dz, d2=params.dx, d3=1.0, d4=params.dt*interval,
           label1=label1, label2=label2, label3=label3, label4=label4,
           title="snapshots", data_format=params.data_format)
end

"""
   append one snapshot to the end of provided stream
"""
function append_one_snapshot(fid::IOStream, spt::Snapshot)

    write(fid, spt.vz); write(fid, spt.vx);
    write(fid, spt.pz); write(fid, spt.px);
    flush(fid)

    return nothing
end

"""
   read one snapshot from binary file in Regularly sampled format
"""
function read_one_snapshot(path::String, it::Int)

    hdr = read_RSheader(path)

    N   = hdr.n1 * hdr.n2
    esize = sizeof(hdr.data_format)

    location = field_location[:data] + (it-1)*N*4*esize
    fid = open(path, "r")
    seek(fid, location)

    vz = zeros(hdr.data_format, N); read!(fid, vz);
    vx = zeros(hdr.data_format, N); read!(fid, vx);
    pz = zeros(hdr.data_format, N); read!(fid, pz);
    px = zeros(hdr.data_format, N); read!(fid, px);
    close(fid)

    return Snapshot(vz, vx, pz, px)
end

"""
   create a header for wavefield in regularly sampled data format
"""
function wavefield_header(params::TdParams, interval::Int64)

    # the size of data
    n1 = params.nz
    n2 = params.nx
    n3 = 3      # include(vz, vx, p)
    n4 = floor(Int64, (params.nt-1)/interval) + 1

    # label
    label1="z-axis"
    label2="x-axis"
    label3="component"
    label4="time step"

    return RegularSampleHeader(n1=n1, n2=n2, n3=n3, n4=n4,
           o1=0.0, o2=0.0, o3=1.0, o4=0.0,
           d1=params.dz, d2=params.dx, d3=1.0, d4=params.dt*interval,
           label1=label1, label2=label2, label3=label3, label4=label4,
           title="wavefield", data_format=params.data_format)
end

"""
   append one wavefield to the end of provided stream
"""
function append_one_wavefield(fid::IOStream, spt::Snapshot, params::TdParams)

    # length of one wavefield
    N = params.nz * params.nx

    # write vz
    for i = 1 : N
        idx = params.spt2wfd[i]
        write(fid, spt.vz[idx])
    end

    # write vx
    for i = 1 : N
        idx = params.spt2wfd[i]
        write(fid, spt.vx[idx])
    end

    # write p (pz+px)
    for i = 1 : N
        idx = params.spt2wfd[i]
        write(fid, spt.pz[idx]+spt.px[idx])
    end
    flush(fid)

    return nothing
end

"""
   save the wavefield at the last time step in a independent file
"""
function write_wavefield(path::String, wfd::Wavefield, params::TdParams)

    # length of wavefield
    n1 = params.nz
    n2 = params.nx
    data_format = eltype(wfd.vz)

    hdr = RegularSampleHeader(n1=n1, n2=n2, n3=3,
                              o1=0.0, o2=0.0, o3=1.0,
                              d1=params.dz, d2=params.dx, d3=1.0,
                              label1="z-axis", label2="x-axis", label3="component",
                              title="wavefield at one time step", data_format=data_format)

    fid = write_RSheader(path, hdr)
    write(fid, wfd.vz)
    write(fid, wfd.vx)
    write(fid, wfd.p )
    close(fid)

    return nothing
end

"""
   read one wavefield at the time step it
"""
function read_one_wavefield(path::String, it::Int64)

    hdr = read_RSheader(path)
    N   = hdr.n1 * hdr.n2
    esize = sizeof(hdr.data_format)

    location = field_location[:data] + (it-1)*N*3*esize
    fid = open(path, "r")
    seek(fid, location)

    vz = zeros(hdr.data_format, N); read!(fid, vz);
    vx = zeros(hdr.data_format, N); read!(fid, vx);
    p  = zeros(hdr.data_format, N); read!(fid, p );
    close(fid)

    return Wavefield(vz, vx, p)
end

"""
   create a header of regularly sampled data for pressure
"""
function pressure_header(params::TdParams, interval::Int64)

    # the size of data
    n1 = params.nz
    n2 = params.nx
    n3 = floor(Int64, (params.nt-1)/interval) + 1

    # label
    label1="z-axis"
    label2="x-axis"
    label3="time step"

    return RegularSampleHeader(n1=n1, n2=n2, n3=n3,
           o1=0.0, o2=0.0, o3=1.0,
           d1=params.dz, d2=params.dx, d3=params.dt*interval,
           label1=label1, label2=label2, label3=label3,
           title="pressure", data_format=params.data_format)
end

"""
   append one pressure to the end of provided stream
"""
function append_one_pressure(fid::IOStream, spt::Snapshot, params::TdParams)

    # length of pressure vector
    N = params.nz * params.nx

    for i = 1 : N
        idx = params.spt2wfd[i]
        tmp = spt.pz[idx] + spt.px[idx]
        write(fid, tmp)
    end
    flush(fid)

    return nothing
end

"""
   read one pressure at a particular time step specified by it
"""
function read_one_pressure(path::String, it::Int64)

    hdr = read_RSheader(path)
    N   = hdr.n1 * hdr.n2
    esize = sizeof(hdr.data_format)

    location = field_location[:data] + (it-1)*N*esize
    fid = open(path, "r")
    seek(fid, location)

    p = zeros(hdr.data_format, N)
    read!(fid, p)

    return reshape(p, hdr.n1, hdr.n2)
end

"""
   Copy the content in spt2 to spt1
"""
function copy_snapshot!(spt1::Snapshot, spt2::Snapshot)

    copyto!(spt1.vz, spt2.vz)
    copyto!(spt1.vx, spt2.vx)
    copyto!(spt1.pz, spt2.pz)
    copyto!(spt1.px, spt2.px)
    return nothing
end

"""
   copy the content in wfd2 to wfd1
"""
function copy_wavefield!(wfd1::Wavefield, wfd2::Wavefield)

    copyto!(wfd1.vz, wfd2.vz)
    copyto!(wfd1.vx, wfd2.vx)
    copyto!(wfd1.p , wfd2.p )
    return nothing
end

"""
   structure save the boundary values of wavefield at all time step. The three
component vz, vx, p are matrix size of nb-by-(nt+1), the one extra time step is used
for born approximation.
"""
struct WavefieldBound{Tv<:AbstractFloat}

    vz :: Union{Matrix{Tv}, Vector{Tv}}
    vx :: Union{Matrix{Tv}, Vector{Tv}}
    p  :: Union{Matrix{Tv}, Vector{Tv}}
end

"""
   Constructor of WavefieldBound, it provide two options, 1: include all time steps, 2: just one time step.
by setting the keyword argument flag=1 or flag=2.
"""
function WavefieldBound(params::TdParams; step_flag=1)

    # length of boundary elements for one time step
    N = length(params.spt2bnd)

    # initialize empty struct
    if step_flag == 1
       return WavefieldBound(zeros(params.data_format, N),
                             zeros(params.data_format, N),
                             zeros(params.data_format, N))

    elseif step_flag == 2
       return WavefieldBound(zeros(params.data_format, N, params.nt),
                             zeros(params.data_format, N, params.nt),
                             zeros(params.data_format, N, params.nt))
    else
       error("non-supported flag, set flag=1 include all time steps, flag=2 just one time step")
    end
end

"""
   Create the  RS header for wavefield boundary
"""
function boundary_header(params::TdParams)

    # the size of data
    n1 = length(params.spt2bnd)
    n2 = 3                      # include(vz, vx, p)

    # the extra boundary value is used to compute time derivative of source-side wavefield
    # via central finite difference
    n3 = params.nt

    # label
    label1="elements"
    label2="component"
    label3="time step"

    return RegularSampleHeader(n1=n1, n2=3, n3=n3,
                               o1=0.0, o2=0.0, o3=1.0, o4=0.0,
                               d1=0.0, d2=0.0, d3=params.dt,
                               label1=label1, label2=label2, label3=label3,
                               title="wavefield boundaries", data_format=params.data_format)
end

"""
   write the boundary values for all time steps
"""
function write_boundary(path_bnd::Ts, bnd::WavefieldBound) where {Ts<:String}

    # create header for wavefield boundaries
    hdr = boundary_header(params)

    # check the size of bnd
    if (hdr.n1, hdr.n3) != size(bnd.vz) || (hdr.n1, hdr.n3) != size(bnd.vx) || (hdr.n1, hdr.n3) != size(bnd.p)
       error("size mismatch")
    end

    # write the header of wavefield boundary
    fid = write_RSheader(path_bnd, hdr)

    # write the data of wavefield boundary
    for i = 1 : hdr.n3
        write(fid, view(bnd.vz, :, i))
        write(fid, view(bnd.vx, :, i))
        write(fid, view(bnd.p , :, i))
    end
    close(fid)

    return nothing
end

"""
   append the wavefield boundary at one time step to the end of fid
"""
function append_one_boundary(fid::IOStream, spt::Snapshot, params::TdParams)

    # length of one wavefield
    N = length(params.spt2bnd)

    # write vz
    for i = 1 : N
        idx = params.spt2bnd[i]
        write(fid, spt.vz[idx])
    end

    # write vx
    for i = 1 : N
        idx = params.spt2bnd[i]
        write(fid, spt.vx[idx])
    end

    # write p (pz+px)
    for i = 1 : N
        idx = params.spt2bnd[i]
        write(fid, spt.pz[idx]+spt.px[idx])
    end
    flush(fid)

    return nothing
end

"""
   read all the boundary value at all time steps
"""
function read_boundary(path::String)

    (hdr, d) = read_RSdata(path)
    return bnd = WavefieldBound(d[:,1,:], d[:,2,:], d[:,3,:])
end

"""
  update the content in WavefieldBound, specifically
"""
function read_one_boundary!(bnd::WavefieldBound, fid::IOStream,
         it::Int64, params::TdParams)

    # byte length of one boundary
    byte_length = length(params.spt2bnd) * sizeof(params.data_format) * 3

    # location of file pointer
    location = field_location[:data] + byte_length * (it-1)
    seek(fid, location)

    # read the boundary value
    read!(fid, bnd.vz)
    read!(fid, bnd.vx)
    read!(fid, bnd.p )

    return nothing
end

"""
   ||snapshot||
"""
function l2norm_snapshot(spt::Snapshot)

    L = zero(eltype(spt.vz))

    L = L + dot(spt.vz, spt.vz)
    L = L + dot(spt.vx, spt.vx)
    L = L + dot(spt.pz, spt.pz)
    L = L + dot(spt.px, spt.px)

    return sqrt(L)
end

"""
   spt = spt1 - spt2
"""
function minus_snapshot(spt1::Snapshot, spt2::Snapshot)

    N = length(spt1.vz)
    data_format = eltype(spt1.vz)

    vz = zeros(data_format, N)
    vx = zeros(data_format, N)
    pz = zeros(data_format, N)
    px = zeros(data_format, N)

    for i = 1 : N
        vz[i] = spt1.vz[i] - spt2.vz[i]
        vx[i] = spt1.vx[i] - spt2.vx[i]
        pz[i] = spt1.pz[i] - spt2.pz[i]
        px[i] = spt1.px[i] - spt2.px[i]
    end

    return Snapshot(vz, vx, pz, px)
end

"""
   spt = spt1 + spt2
"""
function add_snapshot(spt1::Snapshot, spt2::Snapshot)

    N = length(spt1.vz)
    data_format = eltype(spt1.vz)

    vz = zeros(data_format, N)
    vx = zeros(data_format, N)
    pz = zeros(data_format, N)
    px = zeros(data_format, N)

    for i = 1 : N
        vz[i] = spt1.vz[i] + spt2.vz[i]
        vx[i] = spt1.vx[i] + spt2.vx[i]
        pz[i] = spt1.pz[i] + spt2.pz[i]
        px[i] = spt1.px[i] + spt2.px[i]
    end

    return Snapshot(vz, vx, pz, px)
end

"""
   reverse the order of snapshot or wavefield or pressure (path_tmp) and save them in the new
file(path)
"""
function reverse_order(path::String, path_tmp::String; save_flag="snapshot")

    # read header
    hdr = read_RSheader(path_tmp)

    # write header
    fid = write_RSheader(path, hdr)

    # the length of field at one time step
    if save_flag == "snapshot" || save_flag == "wavefield"
       N = hdr.n1 * hdr.n2 * hdr.n3
       nt= hdr.n4
    elseif save_flag == "pressure"
       N = hdr.n1 * hdr.n2
       nt= hdr.n3
    end

    elength = sizeof(hdr.data_format) * N
    fid_tmp = open(path_tmp, "r")
    tmp     = zeros(hdr.data_format, N)

    for it = 1 : nt

        position = field_location[:data] + elength * (nt-it)
        seek(fid_tmp, position)

        read!(fid_tmp, tmp);
        write(fid, tmp)
    end

    close(fid_tmp)
    close(fid)

    return nothing
end
