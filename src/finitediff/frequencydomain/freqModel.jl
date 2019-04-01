type Mesh
    nt    :: Int64
    nz    :: Int64
    nx    :: Int64
    ext   :: Int64
    iflag :: Int64
    dt    :: Float64
    dz    :: Float64
    dx    :: Float64
    function Mesh(nf::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, dz::Float64, dx::Float64)
        return new(nf, nz, nx, ext, iflag, dt, dz, dx)
    end
end

function modelExpand(m::Array{Float64,2}, M::Mesh)
    nz = M.nz; nx = M.nx;
    ext = M.ext; iflag = M.iflag;
    tmpl = repmat(m[:,1]  , 1, ext)
    tmpr = repmat(m[:,end], 1, ext)
    m = hcat(tmpl, m, tmpr)
    if iflag == 1
       tmpu = repmat(m[1,:]'  , ext, 1)
       tmpb = repmat(m[end,:]', ext, 1)
       m = vcat(tmpu, m, tmpb)
    elseif iflag == 2
       tmpb = repmat(m[end,:]', ext, 1)
       m = vcat(m, tmpb)
    end
    return m
end

function getABC(M::Mesh, amp::Float64)
    nz = M.nz; nx = M.nx;
    ext = M.ext; iflag = M.iflag;
    bwd = ((ext:-1:1).^2) ./ ext^2
    fwd = ((1:ext).^2) ./ ext^2
    I2  = (nx+ext+1) : (nx+2*ext)
    Nx  = nx + 2*ext
    if iflag == 1
       Nz = nz + 2*ext
       gamma = zeros(Nz, Nx)
       I1    = (nz+1+ext) : (nz+2*ext)
       gamma[1:ext, :]    += bwd * (ones(Nx))'  # top
       gamma[1:ext,1:ext] -= bwd*bwd'                   # top left corner
       gamma[1:ext,I2]    -= bwd*fwd'                   # top right corner
       gamma[:,1:ext] += ones(Nz) * bwd'              # left
       gamma[:,I2]    += ones(Nz) * fwd'              # right
    elseif iflag == 2
        Nz = nz + ext
        gamma = zeros(Nz, Nx)
        I1    = (nz+1) : (nz+ext)
        gamma[:,1:ext] += ones(nz+ext) * bwd'            # left
        gamma[:,I2]    += ones(nz+ext) * fwd'            # right
    end
    gamma[I1,:]    += fwd * (ones(Nx))'          # bottom
    gamma[I1,1:ext]-= fwd * bwd'                   # bottom left corner
    gamma[I1,I2]-= fwd * fwd'                   # bottom right corner
    gamma *= amp
    return gamma
end

function getSommerfeldBC(M::Mesh, m::Matrix{Float64})
    nz = M.nz ; nx = M.nx;
    ext= M.ext; iflag = M.iflag
    dz = M.dz ; dx = M.dx;
    Nx = nx + 2*ext
    if M.iflag == 1
       Nz = nz + 2*ext
       somm = zeros(Complex128, Nz, Nx)
       somm[1,2:end-1] = -1im * (1.0/dz) .* sqrt(m[1,2:end-1])
    elseif M.iflag == 2
       Nz = nz +   ext
       somm = zeros(Complex128, Nz, Nx)
    end
    somm[end,:] = -1im * (1.0/dz) .* sqrt(m[end,:])
    somm[:  ,1]+= -1im * (1.0/dx) .* sqrt(m[:  ,1])
    somm[:,end]+= -1im * (1.0/dx) .* sqrt(m[:,end])
    return somm
end

function Dxx(n::Int64, h::Float64)
    tmp = 1./(h^2)
    v1 = -tmp * ones(n-1)
    v2 = 2.0 * ones(n)
    v2[1] = 1.
    v2[n] = 1.
    v2    = tmp * v2
    v3    = -tmp * ones(n-1)
    dxx = spdiagm((v1, v2, v3),[-1, 0, 1], n, n)
    return dxx
end


function getLaplacian(M::Mesh)
    nz = M.nz; nx = M.nx;
    dz = M.dz; dx = M.dx;
    ext = M.ext; iflag = M.iflag;
    Nx = nx + 2*ext
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    d2z = Dxx(Nz, dz)
    d2x = Dxx(Nx, dx)
    L   = kron(speye(Nx), d2z)
    L  += kron(d2x, speye(Nz))
    return L
end

function getHelmhotzMtx(L::SparseMatrixCSC{Float64,Int64}, omega::Float64,
                        m::Array{Float64,2}, pml::Array{Float64,2}, somm::Matrix{Complex128})
    size(m) == size(pml) || throw(DimensionMismatch())
    (nz, nx) = size(m)
    mass = zeros(Complex128, nz, nx)
    omega2 = omega * omega
    for i2 = 1 : nx
        for i1 = 1 : nz
            mass[i1,i2] = -omega2 * m[i1,i2] * (1.0 - 1im*pml[i1,i2]/omega) - somm[i1,i2]*omega
        end
    end
    return L + spdiagm(vec(mass))
end

function getAdjointMtx(L::SparseMatrixCSC{Float64,Int64}, omega::Float64,
                        m::Array{Float64,2}, pml::Array{Float64,2}, somm::Matrix{Complex128})
    (n1, n2) = size(m)
    mass = zeros(Complex128, n1, n2)
    omega2 = omega * omega
    for i2 = 1 : n2
        for i1 = 1 : n1
            mass[i1,i2] = -omega2 * m[i1,i2] * (1.0 + 1im*pml[i1,i2]/omega) + somm[i1,i2]*omega
        end
    end
    return L + spdiagm(vec(mass))
end

type FSource
     ns  :: Int64
     isz :: Vector{Int64}
     isx :: Vector{Int64}
     fc  :: Vector{Vector{Complex128}}
     function FSource(isz::Vector{Int64}, isx::Vector{Int64}, f0::Float64, t0::Vector{Float64}, dt::Float64, nf::Int64)
         length(isz) == length(isx) == length(t0) || throw(DimensionMismatch())
         ns = length(isz)
         w = Ricker(f0, dt)
         w = miniPhase(w)
         nw = length(w)
         fc = Vector{Vector{Complex128}}(ns)
         for i = 1 : ns
             tmp = zeros(nf)
             idxl = floor(Int64, t0[i]/dt) + 1
             tmp[idxl:idxl+nw-1] = w
             tmp = fft(tmp)
             fc[i] = tmp
         end
         return new(ns, isz, isx, fc)
     end
end

function miniPhase{T<:Real}(w::Array{T,1})
    nw = length(w)
    nf = 8*nextpow2(nw)
    W = fft(cat(1,w,zeros(nf-nw)))
    A = log(abs(W) + 1.e-8)
    a = 2*ifft(A)
    n2 = floor(Int, nf/2)
    a[n2+2:nf] = 0.0
    a[1] = a[1]/2
    A = exp(fft(a))
    a = real(ifft(A))
    wmin = real(a[1:nw])
end

function src2wfd(src::FSource, iw::Int64, M::Mesh)
    nz = M.nz; nx = M.nx;
    ext= M.ext; iflag = M.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       pz = ext
    elseif iflag == 2
       Nz = nz +   ext
       pz = 0
    end
    Nx = nx + 2*ext;
    q  = zeros(Complex128, Nz, Nx);
    isz = src.isz; isx = src.isx;
    for i = 1 : length(isz)
        q[isz[i]+pz, isx[i]+ext] += src.fc[i][iw]
    end
    return vec(q)
end

"""
get recordings from wavefield
"""
function sample(wfd::Array{Float64,3}, irz::Vector{Int64}, irx::Vector{Int64})
    nr = length(irz)
    (nt, nz, nx) = size(wfd)
    d = zeros(Float64, nt, nr)
    for i = 1 : nr
        d[:,i] = wfd[:,irz[i],irx[i]]
    end
    return d
end
"""
just sample the wavefield of one frequency component
"""
function sample(u::Vector{Complex128}, irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh)
    nz = M.nz; nx = M.nx;
    dz = M.dz; dx = M.dx;
    ext = M.ext; iflag = M.iflag;
    Nx = nx + 2*ext
    if iflag == 1
       Nz = nz + 2*ext
       zu = ext
    elseif iflag == 2
       Nz = nz +   ext
       zu = 0
    end
    u = reshape(u, Nz, Nx)
    nr = length(irz)
    d = zeros(Complex128, nr)
    for i = 1 : nr
        d[i] = u[irz[i]+zu,irx[i]+ext]
    end
    return d
end

"""
inset recordings to wavefield
"""
function insert(d::Vector{Complex128}, irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh)
    nz = M.nz; nx = M.nx;
    dz = M.dz; dx = M.dx;
    ext = M.ext; iflag = M.iflag;
    Nx = nx + 2*ext
    if iflag == 1
       Nz = nz + 2*ext
       zu = ext
    elseif iflag == 2
       Nz = nz +   ext
       zu = 0
    end
    u = zeros(Complex128, Nz, Nx)
    nr = length(irz)
    for i = 1 : nr
        u[irz[i]+zu, irx[i]+ext] = d[i]
    end
    return vec(u)
end

function getDataTime(src::FSource, M::Mesh, pml::Matrix{Float64}, somm::Matrix{Complex128}, L::SparseMatrixCSC{Float64,Int64}, maxiw::Int64; pflag=true)
    wfd = zeros(Complex128, M.nf, M.nz, M.nx)
    df  = 1. / (M.nf * M.dt)
    if M.iflag == 1
       Nz = M.nz + 2*M.ext
       zu = M.ext
    elseif M.iflag == 2
       Nz = M.nz +   M.ext
       zu = 0
    end
    Nx = M.nx + 2*M.ext
    for iw = 2 : maxiw
        if pflag
           println("$iw")
        end
        omega = 2*pi*(iw-1)*df
        H = getHelmhotzMtx(L, omega, m, pml, somm)
        Hfact = lufact(H)
        q = src2wfd(src, iw, M)
        u = Hfact \ q
        u = reshape(u, Nz, Nx)
        tmp = u[zu+1:zu+M.nz, M.ext+1:M.nx+M.ext]
        wfd[iw, :, :] = tmp
        wfd[nf-iw+2, :, :] = conj(tmp)
    end
    return wfd = real(ifft(wfd,1))
end

function intercept(d::Matrix{Float64}, ot::Float64, dt::Float64, nt::Int64)
    il = floor(Int64, ot/dt) + 1
    iu = il + nt - 1
    dcut = d[il:iu, :]
    return dcut
end
