"""
   struct for explosive source
"""
struct Source{Ti<:Int64, Tv<:AbstractFloat}
    isz     :: Ti  # vertical index location
    isx     :: Ti  # horizontal index location
    src2spt :: Ti  # index mapping to snapshot
    src2wfd :: Ti  # index mapping to wavefield
    it_min  :: Ti  # the index lower bound in time axis
    it_max  :: Ti  # the index upper bound in time axis
    p       :: Vector{Tv} # source wavelet
end

"""
   create a ricker wavelet
"""
function ricker(fdoms, dt)

 	  nw = 2.2 / fdom / dt         # decide number of samples
	  nw = 2*floor(Int64, nw/2)+1  # gaurantee the number of samples is odd

	  nc = floor(Int64, nw/2)
	  w  = zeros(nw)               # allocate space for wavelet
	  k  = collect(1:nw)

    for i = 1 : nw
        alpha = (nc-i+1) * fdom * dt * pi
        beta  = alpha^2
        w[i]  = (1 - beta * 2) * exp(-beta)
    end
	  return w
end

"""
  convert wavelet to its corresponding minimum phase version
"""
function convert2miniphase(w::Vector{Tv}) where {Tv<:AbstractFloat}

    nw = length(w)

    # number of samples
    nf = 8*nextpow(2, nw)

    W = fft(vcat(w,zeros(nf-nw)))
    A = log.(abs.(W) .+ 1.e-8)
    a = 2*ifft(A)
    n2 = floor(Int, nf/2)
    a[n2+2:nf] .= 0.0
    a[1] = a[1]/2
    A = exp.(fft(a))
    a = real(ifft(A))
    wmin = real(a[1:nw])

    return wmin
end

"""
   fc is the high cut frequency of low-pass filter, hl is half length of filter
and dt is sampling rate
"""
function tapered_sinc(fc, hl, dt) where {Tv <: AbstractFloat}

    # length of filter
    l = 2*hl
    filter = zeros(l+1)
    ff = fc / (1.0/dt)

    for i = 0 : l
        idx = i - hl
        if idx == 0
           filter[i+1] = 2*pi*ff
        else
           filter[i+1] = sin(2*pi*ff*idx)/idx
        end
        taper = 0.42 - 0.5 * cos(2*pi*i/l) + 0.08*cos(4*pi*i/l)
        filter[i+1] = filter[i+1] * taper
    end

    # normalize filter
    filter = filter / sum(filter)
    return filter
end

"""
   scaled sinc function without tapering, so the decaying of sinc function is slow
"""
function scaled_sinc(fc, hl, dt::Tv) where {Tv<:AbstractFloat}

    l = 2*hl
    filter = zeros(l+1)
    # high cut frequency fraction
    ff = fc / (1.0/dt)

    for i = 0 : l
        idx = i - hl
        if idx == 0
           filter[i+1] = 2*pi*ff
        else
           filter[i+1] = sin(2*pi*ff*idx)/idx
        end
    end

    filter = filter / sum(filter)
    return filter
end

"""
   initialize source structure, we provide four options to specify the source wavelet
   ricker (default option), sinc, miniphase, user provided.
"""
function Source(isz::Ti, isx::Ti, params::ModelParams; ot=0.0, amp=1.0, fdom=20.0, hl=128,
                type_flag="ricker", p=Vector{Float32}(undef,0)) where {Ti<:Int64}

    # non-user provided source wavelet
    if length(p) == 0
       if type_flag == "ricker"
          p = ricker(fdom, params.dt)
          p = amp * p
       elseif flag_type == "miniphase"
          p = ricker(fdom, params.dt)
          p = convert2miniphase(p)
          p = amp * p
       elseif flag_type == "sinc"
          fc = params.data_format(fdom)
          hl = round(Int64, hl)
          p  = tapered_sinc(fc, hl, params.dt)
          p  = amp * p
       end
    end

    # make the data type consistent
    p = convert(Vector{params.data_format}, p)

    # the index time range of source wavelet
    it_min = floor(Int64, ot/params.dt) + 1
    it_max = it_min + nt - 1

    # the auxillary index mapping the location of source to snapshot or wavefield
    src2spt = (isx+params.npml-1)*par.Nz + isz + params.ntop
    src2wfd = (isx-1) * params.nz + isz

    # call the default constructor
    return Source(isz, isx, src2spt, src2wfd, it_min, it_max, p)
end

"""
   get a vector of Source
"""
function get_multi_sources(isz::Vector{Ti}, isx::Vector{Ti}, params::ModelParams; ot=0.0, amp=1.0,
                          fdom=20.0, hl=128, type_flag="ricker", p=Vector{Float32}(undef,0)) where {Ti<:Int64}

    # number of source
    ns = length(isz)
    if length(isx) != ns
       error("length of source coordinate mismatch")
    end

    # starting time
    if typeof(ot) <: Real  # all sources starting at same time
       ot = ot * ones(ns)
    end
    if length(ot) != ns
       error("length of starting time ot mismatch")
    end

    # amplitude for each source
    if typeof(amp) <: Real
       amp = amp * ones(ns)
    end
    if length(amp) != ns
       error("length of amplitude is wrong")
    end

    # dominant frequency for source wavelet
    if typeof(fdom) <: Real
       fdom = fdom * ones(ns)
    end
    if length(fdom) != ns
       error("length of dominant frequency is wrong")
    end

    # sources wavelet(1. p is a Vector{Vector}; 2. p is a vector; 3. p 0 elemenets vector )
    if length(p) == 0  # without user provided wavelet
       if typeof(type_flag) == String               # all sources share same wavelet
          wavelet_type = Vector{String}(undef, ns)
          wavelet_type[:] .= type_flag
          type_flag    = wavelet_type
       elseif typeof(type_flag) == Vector{String}   # use different wavelet
          if length(type_flag) != ns
             error("length of type_flag mismatch")
          end
       end

       # generate the specific source wavelet
       p = Vector{Vector{params.data_format}}(ns)
       for i = 1 : ns
           if type_flag[i] == "ricker"
              tmp = ricker(params.fdom, params.dt)
              p[i]= amp[i] * tmp
           elseif flag_type[i] == "miniphase"
              tmp = ricker(params.fdom, params.dt)
              tmp = convert2miniphase(tmp)
              p[i]= amp[i] * tmp
           elseif flag_type[i] == "sinc"
              f_hc = params.data_format(f_hc)
              hl   = round(Int64, hl)
              tmp  = truncated_sinc(f_hc, hl, params.dt)
              p[i] = amp[i] * tmp
           end
       end

    elseif eltype(p) <: AbstractFloat   # all sources share same given wavelet
       wavelet = Vector{Vector{params.data_format}}(undef, ns)
       for i = 1 : ns
       end

    else                                # sources have different user-give wavelet

    end

    if length(p) != 0
       wavelet = Vector{Vector{params.data_format}}(undef, ns)
       for i = 1 : ns
           wavelet[i] = convert(Vector{params.data_format}, p)
       end

    end

    # allocate memory for vector of sources
    srcs = Vector{Source}(undef, ns)
    for i = 1 : ns
        srcs[i] = Source(isz[i], isx[i], ot[i], amp[i], par; wavelet_type=wavelet_type)
    end
    return srcs
end

"""
   randomly source starting time shift and make the starting time is the integer
times the time interval.
"""
function random_start_time(ns::Ti, dt::Tv, s::Tv) where {Ti<:Int64, Tv<:AbstractFloat}

    ot = s * rand(ns)
    indt = round(Int64, ot/dt)
    ot[:] = indt*dt

    return convert(Vector{typeof(dt)}, ot)
end

"""
   Add source to snapshot(scaled by dt to make them consistent with wave-equation)
"""
function add_source!(spt::Snapshot, src::Source, it::Int64)

    # current times
    tc = (it-1) * src.dt

    # add source
    if src.ot <= tc <= src.tmax

       indt = it - round(Int64, src.ot / src.dt)
       pos  = src.idx1

       if src.vz_flag
          # spt.vz[pos] = spt.vz[pos] + (src.p[indt] * dt) / par.rho[src.isz, src.isx]
          spt.vz[pos] = spt.vz[pos] + (src.p[indt] * src.dt)
       elseif src.vx_flag
          # spt.vx[pos] = spt.vx[pos] + (src.p[indt] * dt) / par.rho[src.isz, src.isx]
          spt.vx[pos] = spt.vx[pos] + (src.p[indt] * src.dt)
       elseif src.p_flag
          # spt.pz[pos] = spt.pz[pos] + src.p[indt] * dt / 2.0
          # spt.px[pos] = spt.px[pos] + src.p[indt] * dt / 2.0
          spt.pz[pos] = spt.pz[pos] + src.p[indt] * src.dt
       else
          error("either add to particle velocity or stress")
       end
    end
    return nothing
end

"""
   add source to wavefield
"""
function add_source!(wfd::Wavefield, src::Source, it::Int64)

    tc = (it-1) * src.dt

    if src.ot <= tc <= src.tmax
       indt = it - round(Int64, src.ot / src.dt)
       pos  = src.idx2
       wfd.p[pos] = wfd.p[pos] + src.p[indt] * src.dt
    end
    return nothing
end

"""
   subtract source from wavefield
"""
function subtract_source!(wfd::Wavefield, src::Source, it::Int64)

    tc = (it-1) * src.dt

    if src.ot <= tc <= src.tmax
       indt = it - round(Int64, src.ot / src.dt)
       pos  = src.idx2
       wfd.p[pos] = wfd.p[pos] - src.p[indt] * src.dt
    end
    return nothing
end

# inject multiple sources simultaneously
function add_multi_sources!(spt::Snapshot, srcs::Vector{Source}, it::Int64)

    # number of sources
    ns = length(srcs)

    for i = 1 : ns
        add_source!(spt, srcs[i], it)
    end
    return nothing
end

"""
   find the time range of multiple sources
"""
function time_range_multisources(srcs::Vector{Source})

    tmax = 0.0
    tmin = 0.0
    ns   = length(srcs)

    for i = 1 : ns
        tmp1 = srcs[i].ot
        tmp2 = srcs[i].tmax

        if tmp1 < tmin
           tmin = tmp1
        end
        if tmp2 > tmax
           tmax = tmp2
        end
    end

    return tmin, tmax
end



# ========================need to be updated ===================================
# """
#    generate a snapshot via multiple sources
# """
# function srcs2spt(srcs::Array{Source,1}, it::Int64)
#
#     nz = srcs[1].nz ; nx = srcs[1].nx;
#     ext= srcs[1].ext; iflag=srcs[1].iflag;
#     dt = srcs[1].dt
#     spt = InitSnapShot(nz, nx, ext, iflag, dt, it)
#
#     AddMultiSources!(spt, srcs)
#     return spt
# end
#
# """
#    generate wavefield via multiple sources and write the snapshots to hard drive
# """
# function WriteSrcs2SnapShots(path::String , srcs::Array{Source,1})
#     nz = srcs[1].nz
#     (tl, tu) = SrcRange(srcs)
#     nt  = round(Int64, (tu-tl)/dt) + 1
#     fid = open(path, "w")
#     write(fid, nz, nx, ext, iflag, dt)
#     for it = 1 : nt
#         it1 = round(Int64, (tl+(it-1)*dt)/dt+1)
#         spt = srcs2spt(srcs, it1)
#         writeSnapShot(fid, spt)
#     end
#     close(fid)
#     return nothing
# end
#
# """
#    extract a wavelet from the stress field
# """
# function multiStress2wlet(path::String , iz::Int64, ix::Int64)
#     fid = open(path, "r")
#     (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
#     w = zeros(nt)
#     for it = 1 : nt
#         sts= readStress(path, it)
#         w[it] = sts.p[iz,ix]
#     end
#     return w
# end
