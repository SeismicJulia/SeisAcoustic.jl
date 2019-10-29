function trace_interpolation(d::Array{Tv}, order::Ti) where {Tv<:AbstractFloat, Ti<:Int64}

    Ndims = ndims(d)
    Ndims == 1 ? nx = 1 :  nx = prod(size(d)[2:Ndims])
    dims  = size(d)
    N     = dims[1]
    Npad  = order*N
    nf    = 2*Npad
    nw    = convert(Int,floor(nf/2)) + 1
    di    = zeros(eltype(d),nf,dims[2:end]...)

    for k = 1:nx

        dd = zeros(eltype(d),nf)
        dd[1:order:Npad] = d[1:1:end,k]
        D = fft(dd)
        nyq = convert(Int,floor(nw/order)) + 1
        for iw = nyq : nw
            D[iw] *= 0
        end

        # symmetries
        for iw=nw+1:nf
            D[iw] = conj(D[nf-iw+2])
        end
        di[:,k] = real(ifft(D,1))
    end

    return  reshape(di[1:Npad,:]*order,Npad,dims[2:end]...)
end
