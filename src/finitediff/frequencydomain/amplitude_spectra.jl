"""
   flag = true, take the traditional convention for Fourier transform
   flag = false, take the opposite convertion
   fhigh_end = 0.0, the high end of frequency axis is Nyquensit frequency
"""
function amplitude_spectra(f::Vector{Tv}, dt::Tv, nf::Ti;
         flag=true, fhigh_end=0.0) where {Tv<:AbstractFloat, Ti<:Int64}

    nt   = length(f)

    # padding zeros to input
    if nf < nt
       nf = nt
    else
       fp = vcat(f, zeros(typeof(dt), nf-nt))
    end

    #  make sure nf is an even number
    if mod(nf, 2) != 0
       nf = nf + 1
       fp = vcat(fp, zero(typeof(dt)))
    end

    fmax = 1.0 / dt
    df = fmax / nf

    # Fourier transform, flag is a indicator for convention
    # taking normalized transform
    if flag
       F = 1.0 / sqrt(nf) * fft(fp)
    else
       F = sqrt(nf) * ifft(fp)
    end

    # only output amplitude spectra
    if fhigh_end == 0.0
       nw = floor(Int64, nf/2) + 1
    else
       nw = floor(Int64, fhigh_end/df) + 1
    end

    faxis = collect(0:nw) * df
    amp= abs.(F[1:nw+1])

    return faxis, amp

end


function amplitude_spectra(d::Matrix{Tv}, dt::Tv, nf::Ti;
         flag=true, fhigh_end=0.0) where {Tv<:AbstractFloat, Ti<:Int64}

    (nt, ntrace) = size(d)

    # padding zeros to input
    if nf < nt
       nf = nt
    else
       dp = vcat(d, zeros(typeof(dt), nf-nt, ntrace))
    end

    #  make sure nf is an even number
    if mod(nf, 2) != 0
       nf = nf + 1
       dp = vcat(dp, zero(typeof(dt), 1, ntraces))
    end

    fmax = 1.0 / dt
    df = fmax / nf

    # Fourier transform, flag is a indicator for convention
    # taking normalized transform
    if flag
       F = 1.0 / sqrt(nf) * fft(dp, 1)
    else
       F = sqrt(nf) * ifft(fp)
    end

    # only output amplitude spectra
    if fhigh_end == 0.0
       nw = floor(Int64, nf/2) + 1
    else
       nw = floor(Int64, fhigh_end/df) + 1
    end
    faxis = collect(0:nw) * df

    # smooth the amplitude spectra by stacking
    amp = abs.(F[1:nw+1, :])
    amp = sum(amp, 2)

    return faxis, amp

end


function time_derivative(f::Vector{Tv}, dt::Tv; nf=1024) where {Tv<:AbstractFloat}

    nt = length(f)

    # padding zeros to input
    if nf < nt
       nf = nt
    else
       fp = vcat(f, zeros(typeof(dt), nf-nt))
    end

    #  make sure nf is an even number
    if mod(nf, 2) != 0
       nf = nf + 1
       fp = vcat(fp, zero(typeof(dt)))
    end

    fmax = 1.0 / dt
    df   = fmax / nf

    F = 1.0 / sqrt(nf) * fft(fp)
    F[1] = 0.0+0.0im
    for i = 2 : floor(Int64, nf/2)
        omega = (i-1) * 2 * pi * df
        F[i]  = 1im * omega * F[i]
        F[nf-i+2] = conj(F[i])
    end
    return sqrt(nf) * real(ifft(F))[1:nt]
end

"""
   band pass filtering along the first dimension
"""
function band_pass_filter(din::Array{Tv}, dt::Tv,
         f1::Tv, f2::Tv, f3::Tv, f4::Tv) where {Tv<:AbstractFloat}

    # get the data format
    data_format = eltype(din)

    # dimensions of input data
    dims = collect(size(din))

    # padding zeros
    nt  = dims[1]
    nf  = 4 * nextpow2(nt)
    dims_pad = copy(dims)
    dims_pad[1] = nf - nt
    pad = zeros(data_format, dims_pad...)
    dp  = vcat(din, pad)
    dp  = fft(dp, 1)

    # the index of Nyquensit frequency
    nyqind = floor(Int64, nf/2)+1

    # the index of corner frequency
    i1 = floor(Int64, f1*dt*nf) + 1
    i2 = floor(Int64, f2*dt*nf) + 1
    i3 = floor(Int64, f3*dt*nf) + 1
    i4 = floor(Int64, f4*dt*nf) + 1

    # amplitude taper in frequency domain
    up  = collect(1:(i2-i1)) / (i2-i1)
    down= collect((i4-i3):-1:1) / (i4-i3)
    aux = vcat(zeros(i1), up, ones(i3-i2), down, zeros(nyqind-i4))
    aux2= flipdim(aux[2:nyqind-1], 1)

    # complex band pass filter
    c   = 0.
    amp = vcat(aux, aux2);
    pha = (pi/180.) * vcat(0., c*ones(nyqind-2), 0., -c*ones(nyqind-2))
    bpfilter = amp .* exp.(im*pha)

    # applying filter trace by trace
    if length(dims) == 1
       dout  = zeros(Complex{data_format}, nf)
       dout  = bpfilter .* dp
       dout  = ifft(dout, 1)
       dfilt = real(dout[1:nt])

    elseif length(dims) >= 2
       nx   = prod(dims[2:end])
       dout = zeros(Complex{data_format}, nf, nx)
       dp   = reshape(dp, nf, nx)
       for ix = 1: nx
           dout[:, ix] = bpfilter .* dp[:, ix]
       end
       dout  = ifft(dout, 1)
       dfilt = real(dout[1:nt, :])
       dfilt = reshape(dfilt, dims...)
    end

    return dfilt
end

# upsampling of time series via sinc interpolation
function sinc_resample1D(d::Vector{Tv}, order::Ti) where {Tv<:AbstractFloat, Ti<:Integer}

    N    = length(d)
	  Npad = order*N - 1
    nf   = 2*Npad
    nw   = convert(Int,floor(nf/2)) + 1
	  dd   = zeros(eltype(d), nf)
    dd[1:order:Npad] = d[1:1:end]
	  D = fft(dd)
    nyq = convert(Int,floor(nw/order)) + 1

	  for iw = nyq : nw
		    D[iw] *= 0
	  end

	  # symmetries
	  for iw=nw+1:nf
		    D[iw] = conj(D[nf-iw+2])
	  end

    dd = real(ifft(D,1))
	  return dd[1:Npad]*order

end
