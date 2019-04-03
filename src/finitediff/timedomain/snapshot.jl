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
   Initialize a empty snapshot, all the fields are zero vectors
"""
function Snapshot(params::ModelParams)

    data_format = params.data_format
    N           = params.Nz * params.Nx

    return Snapshot(zeros(data_format,N), zeros(data_format,N),
                    zeros(data_format,N), zeros(data_format,N))
end

"""
   Struct for wavefiled which include three fields(vz, vx, p), PML boundary part are cropped
"""
struct Wavefield{Tv<:AbstractFloat}
     vz :: Vector{Tv}
     vx :: Vector{Tv}
     p  :: Vector{Tv}
end

"""
   Initialize zero wavefield
"""
function Wavefield(params::ModelParams)

    data_format = params.data_format
    N           = params.nz * params.nx

    return Wavefield(zeros(data_format,N), zeros(data_format,N), zeros(data_format,N))
end

"""
   Crop the PML layers of snapshot and sum pz, px to obtain wavefield at one time step
"""
function sample_spt2wfd(spt::Snapshot, params::ModelParams)

    # initial a empty wavefield
    wfd = Wavefield(params)

    i = 0
    for j in params.spt2wfd
        i = i + 1
        wfd.vz[i] = spt.vz[j]
        wfd.vx[i] = spt.vx[j]
        wfd.p[i]  = spt.pz[j] + spt.px[j]
    end
    return wfd
end

"""
   Create header of regularly sampled data for snapshot
"""
function snapshot_header(params::ModelParams)

    # the size of data
    n1 = params.Nz
    n2 = params.Nx
    n3 = 4       # include vz, vx, pz, px
    n4 = params.nt

    # label
    label1="z-axis"
    label2="x-axis"
    label3="component"
    label4="time step"

    return RegularSampleHeader(n1=n1, n2=n2, n3=n3, n4=n4,
           o1=0.0      , o2=0.0      , o3=1.0, o4=0.0      ,
           d1=params.dz, d2=params.dx, d3=1.0, d4=params.dt,
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
function wavefield_header(params::ModelParams)

    # the size of data
    n1 = params.nz
    n2 = params.nx
    n3 = 3      # include(vz, vx, p)
    n4 = params.nt

    # label
    label1="z-axis"
    label2="x-axis"
    label3="component"
    label4="time step"

    return RegularSampleHeader(n1=n1, n2=n2, n3=n3, n4=n4,
           o1=0.0, o2=0.0, o3=1.0, o4=0.0,
           d1=params.dz, d2=params.dx, d3=1.0, d4=params.dt,
           label1=label1, label2=label2, label3=label3, label4=label4,
           title="wavefield", data_format=params.data_format)
end

"""
   append one wavefield to the end of provided stream
"""
function append_one_wavefield(fid::IOStream, spt::Snapshot, params::ModelParams)

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
function write_wavefield(path::String, wfd::Wavefield, params::ModelParams)

    # length of wavefield
    n1 = params.nz
    n2 = params.nx
    data_format = eltype(wfd.vz)

    hdr = RegularSampleHeader(n1=n1, n2=n2, n3=3,
                              o1=1 , o2=1 , o3=1,
                              d1=1 , d2=1 , d3=1,
                              label1="wavefield value", label2="component",
                              title="last wavefield", data_format=data_format)

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
function pressure_header(params::ModelParams)

    # the size of data
    n1 = params.nz
    n2 = params.nx
    n3 = params.nt

    # label
    label1="z-axis"
    label2="x-axis"
    label3="time step"

    return RegularSampleHeader(n1=n1, n2=n2, n3=n3,
           o1=0.0, o2=0.0, o3=1.0,
           d1=params.dz, d2=params.dx, d3=params.dt,
           label1=label1, label2=label2, label3=label3,
           title="pressure", data_format=params.data_format)
end

"""
   append one pressure to the end of provided stream
"""
function append_one_pressure(fid::IOStream, spt::Snapshot, params::ModelParams)

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
function WavefieldBound(params::ModelParams; step_flag=1)

    # length of boundary elements for one time step
    N = length(params.spt2bnd)

    # the one extra time step is used for computing the time derivative of the source-side wavefield
    nt= params.nt + 1

    # initialize empty struct
    if flag == 1
       return WavefieldBound(zeros(params.data_format, N),
                             zeros(params.data_format, N),
                             zeros(params.data_format, N))

    elseif flag == 2
       return WavefieldBound(zeros(params.data_format, N, nt),
                             zeros(params.data_format, N, nt),
                             zeros(params.data_format, N, nt))
    else
       error("non-supported flag, set flag=1 include all time steps, flag=2 just one time step")
    end
end

"""
   Create the  RS header for wavefield boundary
"""
function boundary_header(params::ModelParams)

    # the size of data
    n1 = length(params.spt2bnd)
    n2 = 3                      # include(vz, vx, p)

    # the extra boundary value is used to compute time derivative of source-side wavefield
    # via central finite difference
    n3 = params.nt+1

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
function append_one_boundary(fid::IOStream, spt::Snapshot, params::ModelParams)

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
function read_one_boundary!(bnd::WavefieldBound, fid::IOStream)

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

# """
#    reverse the order of snapshot or wavefield or pressure (path_tmp) and save them in the new
# file(path)
# """
# function reverse_order(path::String, path_tmp::String; data_type="snapshot")
#
#     # read header
#     hdr = read_RSheader(path_tmp)
#
#     # write header
#     fid = write_RSheader(path, hdr)
#
#     # the length of field at one time step
#     if data_type == "snapshot" || data_type == "wavefield"
#        N = hdr.n1 * hdr.n2 * hdr.n3
#     elseif datatype == "pressure"
#        N = hdr.n1 * hdr.n2
#     end
#
#     elength = sizeof(hdr.data_format) * N
#     fid_tmp = open(path_tmp, "r")
#     tmp     = zeros(hdr.data_format, N)
#
#     for it = 1 : hdr.nt
#
#         position = field_location[:data] + elength * (hdr.nt-it)
#         seek(fid_tmp, position)
#
#         read!(fid_tmp, tmp);
#         write(fid, tmp)
#     end
#
#     close(fid_tmp)
#     close(fid)
#
#     return nothing
# end

#////////////////////// fixed them at this step ////////////////////////////////
# """
#    compute the l1l2 mixed norm (group sparsity)
# """
# function l1l2norm_pressure(path::String)
#
#     hdr = read_USdata(path_tmp, data_flag=false)
#     tmp =
#
#
#     for it = 1 : hdr.n2
#         sts = read(fid, Float64, nz*nx)
#         tmp = tmp + sts.^2
#     end
#     close(fid)
#     tmp = sqrt(tmp)
#     return tmp
# end

# function L2normStress(path::String )
#     (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
#     fid = open(path, "r")
#     seek(fid, sizeof(Float64)*5)
#     tmp = 0.0
#     for it = 1 : nt
#         spt = read(fid, Float64, nz*nx)
#         tmp = tmp + dot(spt, spt)
#     end
#     close(fid)
#     tmp = sqrt(tmp)
#     return tmp
# end
#
# function scaleStress(path::String , alpha::Float64)
#     (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
#     fid = open(path, "r+")
#     pre = sizeof(Float64) * 5
#     lspt = nz*nx*sizeof(Float64)
#     for it = 1 : nt
#         position = pre + (it-1)*lspt
#         seek(fid, position)
#         spt = read(fid, Float64, nz*nx)
#         spt[:] = spt * alpha
#         seek(fid, position)
#         write(fid, spt)
#     end
#     close(fid)
#     return nothing
# end
#
# function addStress2Spt!(spt::SnapShot, path::String )
#     it = spt.it
#     nz = spt.nz;   nx = spt.nx;
#     ext = spt.ext; iflag = spt.iflag;
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nz +   ext
#     end
#     Nx = nx + 2*ext
#     sts = readStress(path, it)
#     p   = sts.p * 1/2
#     tmp = reshape(spt.pz, Nz, Nx)
#     tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + p
#     spt.pz = vec(tmp)
#     tmp = reshape(spt.px, Nz, Nx)
#     tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + p
#     spt.px = vec(tmp)
#     return nothing
# end
#
# function addReflection2Spt!(spt::SnapShot, I::Array{Float64,2}, path::String )
#     it = spt.it
#     nz = spt.nz; nx = spt.nx;
#     ext = spt.ext; iflag = spt.iflag;
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nz +   ext
#     end
#     Nx = nx + 2*ext
#     sts = readStress(path, it)
#     ref = sts.p .* I * 1/2
#     tmp = reshape(spt.pz, Nz, Nx)
#     tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + ref
#     spt.pz[:] = vec(tmp)
#     tmp = reshape(spt.px, Nz, Nx)
#     tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + ref
#     spt.px[:] = vec(tmp)
#     return nothing
# end
#
# function spt2dis(w::Array{Float64,1}, spt::SnapShot)
#     nz = spt.nz ; nx = spt.nx;
#     ext= spt.ext; iflag = spt.iflag;
#     it = spt.it
#     if length(w) < spt.it
#        error("out of source time range")
#     end
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nz +   ext
#     end
#     Nx = nx + 2*ext
#     pz = reshape(spt.pz, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
#     px = reshape(spt.px, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
#     dis = 1/2 * w[it] * (pz+px)
#     return dis
# end
#
# function spt2wlet(w::Array{Float64,1}, dis::Array{Float64,2}, spt::SnapShot)
#     nz = spt.nz ; nx = spt.nx;
#     ext= spt.ext; iflag = spt.iflag;
#     it = spt.it
#     if length(w) < it
#        error("out of source time range")
#     end
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nz +   ext
#     end
#     Nx = nx + 2*ext
#     p = reshape(spt.pz, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx] + reshape(spt.px, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
#     p = vec(p)
#     w[it] = 1/2 * dot(vec(dis), vec(p))
#     return nothing
# end
#
# function diswlet2spt(dis::Array{Float64,1}, w::Array{Float64,1}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
#     if length(w) < it
#        error("out of source time range")
#     end
#     if length(dis) != nz*nx
#        error("length of dis does not much model size")
#     end
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nz +   ext
#     end
#     Nx = nx + 2*ext
#     spt = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
#     tmp = zeros(Nz, Nx); dis_tmp = reshape(dis, nz, nx)
#     tmp[zupper+1: zupper+nz, ext+1:ext+nx] = 1/2 * w[it] * dis_tmp
#     tmp = vec(tmp)
#     spt.pz[:] = tmp[:]
#     spt.px[:] = tmp[:]
#     return spt
# end
