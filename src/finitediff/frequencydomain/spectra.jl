"""
   compute the frequency amplitude spectra for a vector signal, if fhigh_end is not
specified, the default is nyquensit frequency.
"""
function amplitude_spectra(f::Vector{Tv}, dt;
                           nf=0, fhigh_end=0.0) where {Tv<:AbstractFloat, Ti<:Int64}

    # length of signal
    nt = length(f)

    # determine the number of zeros need to be padded
    if nf == 0 || nf < nt
       nf = nextpow(2, nt)
    end

    #  make sure nf is an even number
    if mod(nf, 2) != 0
       nf = nf + 1
    end

    # padding zeros
    f_pad = vcat(f, zeros(Tv, nf-nt))

    # frequency sampling interval
    fmax = 1.0 / dt
    df   = fmax/ nf

    # taking normalized transform
    F = 1.0 / sqrt(nf) * fft(f_pad)

    # only output amplitude spectra
    if fhigh_end == 0.0
       nw = floor(Int64, nf/2) + 1
    else
       nw = floor(Int64, fhigh_end/df) + 1
    end

    # make the frequency axis
    faxis = collect(0 : nw-1) * df

    # make amplitude spectrum
    amp= abs.(F[1:nw])

    return faxis, amp

end

"""
   compute the mean time-frequency spectrum when input is a matrix, the amplitude spectra
of each column is averaged.
"""
function amplitude_spectra(d::Matrix{Tv}, dt;
         nf=0, fhigh_end=0.0) where {Tv<:AbstractFloat, Ti<:Int64}

    # size of input
    (nt, ntrace) = size(d)

    # determine the number of zeros need to be padded
    if nf == 0 || nf < nt
       nf = nextpow(2, nt)
    end

    #  make sure nf is an even number
    if mod(nf, 2) != 0
       nf = nf + 1
    end

    # padding zeros to input
    d_pad = vcat(d, zeros(Tv, nf-nt, ntrace))

    # frequency sampling interval
    fmax = 1.0 / dt
    df   = fmax/ nf

    # taking normalized transform
    F = 1.0 / sqrt(nf) * fft(d_pad, 1)

    # only output amplitude spectra
    if fhigh_end == 0.0
       nw = floor(Int64, nf/2) + 1
    else
       nw = floor(Int64, fhigh_end/df) + 1
    end
    faxis = collect(0:nw-1) * df

    # smooth the amplitude spectra
    amp = abs.(F[1:nw, :])
    amp = sum(amp; dims=2) / ntrace

    return faxis, amp
end

# """
#    compute time derivative in frequency domain
# """
# function spectra_time_derivative(f::Vector{Tv}, dt) where {Tv<:AbstractFloat}
#
#     nt = length(f)
#
#     # padding zeros to input
#     nf = nextpow(2, nt)   # nf is always an even number
#     fp = vcat(f, zeros(Tv, nf-nt))
#
#     # frequency sampling interval
#     fmax = 1.0 / dt
#     df   = fmax/ nf
#
#     # transform to frequency domain
#     F = fft(fp)
#
#     # index of nyquensit frequency
#     iw_upper = floor(Int64, nf/2)
#
#     # multiplication
#     F[1] = zero(Complex{Tv})
#     for i = 2 : iw_upper
#         omega = 2 * pi * (i-1) * df
#         F[i]  = 1im * omega * F[i]
#
#         # conjugcy of real signal
#         F[nf-i+2] = conj(F[i])
#     end
#
#     # inverse fft, take the real part and truncate to original length
#     return real(ifft(F))[1:nt]
#
# end
#
# """
#    band pass filtering along the first dimension
# """
# function band_pass_filter(din::Array{Tv}, dt::Tv,
#          f1::Tv, f2::Tv, f3::Tv, f4::Tv) where {Tv<:AbstractFloat}
#
#     # get the data format
#     data_format = eltype(din)
#
#     # dimensions of input data
#     dims = collect(size(din))
#
#     # padding zeros
#     nt  = dims[1]
#     nf  = 4 * nextpow2(nt)
#     dims_pad = copy(dims)
#     dims_pad[1] = nf - nt
#     pad = zeros(data_format, dims_pad...)
#     dp  = vcat(din, pad)
#     dp  = fft(dp, 1)
#
#     # the index of Nyquensit frequency
#     nyqind = floor(Int64, nf/2)+1
#
#     # the index of corner frequency
#     i1 = floor(Int64, f1*dt*nf) + 1
#     i2 = floor(Int64, f2*dt*nf) + 1
#     i3 = floor(Int64, f3*dt*nf) + 1
#     i4 = floor(Int64, f4*dt*nf) + 1
#
#     # amplitude taper in frequency domain
#     up  = collect(1:(i2-i1)) / (i2-i1)
#     down= collect((i4-i3):-1:1) / (i4-i3)
#     aux = vcat(zeros(i1), up, ones(i3-i2), down, zeros(nyqind-i4))
#     aux2= flipdim(aux[2:nyqind-1], 1)
#
#     # complex band pass filter
#     c   = 0.
#     amp = vcat(aux, aux2);
#     pha = (pi/180.) * vcat(0., c*ones(nyqind-2), 0., -c*ones(nyqind-2))
#     bpfilter = amp .* exp.(im*pha)
#
#     # applying filter trace by trace
#     if length(dims) == 1
#        dout  = zeros(Complex{data_format}, nf)
#        dout  = bpfilter .* dp
#        dout  = ifft(dout, 1)
#        dfilt = real(dout[1:nt])
#
#     elseif length(dims) >= 2
#        nx   = prod(dims[2:end])
#        dout = zeros(Complex{data_format}, nf, nx)
#        dp   = reshape(dp, nf, nx)
#        for ix = 1: nx
#            dout[:, ix] = bpfilter .* dp[:, ix]
#        end
#        dout  = ifft(dout, 1)
#        dfilt = real(dout[1:nt, :])
#        dfilt = reshape(dfilt, dims...)
#     end
#
#     return dfilt
# end
#
# # upsampling of time series via sinc interpolation
# function sinc_resample1D(d::Vector{Tv}, order::Ti) where {Tv<:AbstractFloat, Ti<:Integer}
#
#     N    = length(d)
# 	  Npad = order*N - 1
#     nf   = 2*Npad
#     nw   = convert(Int,floor(nf/2)) + 1
# 	  dd   = zeros(eltype(d), nf)
#     dd[1:order:Npad] = d[1:1:end]
# 	  D = fft(dd)
#     nyq = convert(Int,floor(nw/order)) + 1
#
# 	  for iw = nyq : nw
# 		    D[iw] *= 0
# 	  end
#
# 	  # symmetries
# 	  for iw=nw+1:nf
# 		    D[iw] = conj(D[nf-iw+2])
# 	  end
#
#     dd = real(ifft(D,1))
# 	  return dd[1:Npad]*order
#
# end
