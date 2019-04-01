function window(woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                M::Mesh, u::Vector{Complex128}; iflag="fwd")
    nz = M.nz; nx = M.nx; ext = M.ext;
    if M.iflag == 1
       Nz = nz + 2*ext
       zu = ext
    elseif M.iflag == 2
       Nz = nz +   ext
       zu = 0
    end
    Nx = nx + 2*ext
    if iflag == "fwd"
       q = zeros(Complex128, Nz, Nx)
       izl = zu+woz ; izu = izl + wnz - 1;
       ixl = ext+wox; ixu = ixl + wnx - 1;
       q[izl:izu, ixl:ixu] = reshape(u, wnz, wnx);
    elseif iflag == "adj"
       izl = zu+woz ; izu = izl + wnz - 1;
       ixl = ext+wox; ixu = ixl + wnx - 1;
       q = reshape(u, Nz, Nx)[izl:izu, ixl:ixu]
    end
    return vec(q)
end

function forwardOperator(s::Vector{Complex128},
                         woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                         irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh,
                         H::Base.SparseArrays.UMFPACK.UmfpackLU{Complex128,Int64})
    q = window(woz, wox, wnz, wnx, M, s, iflag="fwd")
    u = H \ q
    d = sample(u, irz, irx, M)
    return d
end

function adjointOperator(d::Vector{Complex128},
                          woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                          irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh,
                          A::Base.SparseArrays.UMFPACK.UmfpackLU{Complex128,Int64})
    u = insert(d, irz, irx, M)
    q = A \ u
    s = window(woz, wox, wnz, wnx, M, q, iflag="adj")
    return s
end

function adjoint_stack(d::Matrix{Float64}, dt::Float64, miniw::Int64, maxiw::Int64,
                  L::SparseMatrixCSC{Float64,Int64}, m::Matrix{Float64},
                  pml::Matrix{Float64}, somm::Matrix{Complex128},
                  woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                  irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh)
    fd = fft(d,1)
    (nf, nr) = size(fd)
    df = 1/(dt*nf)
    loc= zeros(wnz*wnx)
    for iw = miniw : maxiw
         omega = 2*pi*(iw-1)*df
         dobs  = fd[iw,:]
         H = getAdjointMtx(L, omega, m, pml, somm)
         H = lufact(H)
         tmp = adjointOperator(dobs, woz, wox, wnz, wnx, irz, irx, M, H)
         loc = loc .+ abs(tmp)
         println("$iw")
    end
    return reshape(loc, wnz, wnx)
end

function CG_stack(d::Matrix{Float64}, dt::Float64, maxiw::Int64,
                  L::SparseMatrixCSC{Float64,Int64}, m::Matrix{Float64}, pml::Matrix{Float64},
                  woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                  irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh;
                  maxiter=5, iw_increment=6)
    fd = fft(d,1)
    (nf, nr) = size(fd)
    df = 1/(dt*nf)
    loc= zeros(wnz*wnx)
    for iw = 2 : iw_increment : maxiw
         omega = 2*pi*(iw-1)*df
         dobs  = fd[iw,:]
         H = getHelmhotzMtx(L, omega, m, pml); H = lufact(H);
         A = getAdjointMtx(L, omega, m, pml);  A = lufact(A);
         (tmp, J) = CG(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A);
         loc = loc .+ (abs(tmp)).^2
         println("$iw")
    end
    return reshape(loc, wnz, wnx)
end


"""
power's method to determine the maximum eigenvalue of H * H'(H'=A)
"""
function fpower(woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                irz::Vector{Int64},irx::Vector{Int64}, M::Mesh,
                H::Base.SparseArrays.UMFPACK.UmfpackLU{Complex128,Int64},
                A::Base.SparseArrays.UMFPACK.UmfpackLU{Complex128,Int64};
                maxiter=60, pflag=false)
    q = rand(Complex128, wnz*wnx)
    lambda = 0.
    for i = 1 : maxiter
        d = forwardOperator(q, woz, wox, wnz, wnx, irz, irx,  M, H)
        y = adjointOperator(d, woz, wox, wnz, wnx, irz, irx,  M, A)
        v = norm(q)
        q = y / v
        lambda = v
        if pflag == true
           println("iteration: $i, lambda=$lambda")
        end
    end
    return lambda
end

"""
fista method for source location in frequency domain
"""
function ffista(d::Vector{Complex128}, woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh,
                H::Base.SparseArrays.UMFPACK.UmfpackLU{Complex128, Int64},
                A::Base.SparseArrays.UMFPACK.UmfpackLU{Complex128, Int64},
                mu::Float64, lambda::Float64; maxiter=100, pflag=false)
    x = zeros(Complex128, wnz*wnx)
    T = mu / (2*lambda)
    t = 1
    yk = copy(x)
    J  = zeros(maxiter)
    for k = 1 : maxiter
        tmpx = copy(x)
        dcal = forwardOperator(yk, woz, wox, wnz, wnx, irz, irx, M, H)
        r = dcal - d;
        cost = norm(r) + mu * sum(abs(yk)); J[k] = cost;
        if pflag == true
           println("iteration $k, cost $cost")
        end
        x = adjointOperator(r, woz, wox, wnz, wnx, irz, irx, M, A)
        x = yk - x/lambda;
        fsoftThresh!(x, T);
        tmpt = t
        t = (1 + sqrt(1+4*t^2)) / 2
        yk = x + (tmpt-1)/t * (x-tmpx)
    end
    return x, J
end

function fsoftThresh!(x::Vector{Complex128}, t::Float64)
    tmp = abs(x) - t
    for i = 1 : length(x)
        z = abs(x[i])
        u = maximum([z-t, 0.]);
        l = z + 1.e-16
        x[i] = (u/l) * x[i]
    end
    return nothing
end

function CG(d::Vector{Complex128}, woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
            irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh,
            H::Base.SparseArrays.UMFPACK.UmfpackLU{Complex{Float64},Int64},
            A::Base.SparseArrays.UMFPACK.UmfpackLU{Complex{Float64},Int64};
            Niter=10, mu=0, tol=1.0e-10)

    cost = Float64[]
    r = copy(d)
    g = adjointOperator(r, irz, irx, woz, wox, wnz, wnx, M, A)
    m = zeros(g)
    s = copy(g)
    gamma = (norm(g))^2
    gamma0= gamma
    cost0 = (norm(r))^2
    push!(cost,1.0)
    for iter = 1 : Niter
        t = forwardOperator(s, irz, irx, woz, wox, wnz, wnx, M, H)
	      delta = (norm(t))^2 + mu * (norm(s))^2
	      alpha = gamma/delta
	      m = m + alpha * s
	      r = r - alpha * t
	      g = adjointOperator(r, irz, irx, woz, wox, wnz, wnx, M, A)
	      g = g - mu * m
	      gamma0 = gamma
	      gamma  = (norm(g))^2
        cost1  = (norm(r))^2 + mu * (norm(m))^2
        tmp    = cost1 / cost0
        println("iteration $iter, cost $tmp")
        push!(cost,tmp)
	      beta = gamma / gamma0
	      s    = beta * s + g
	      # if (sqrt(gamma) <= sqrt(gamma00) * tol)
	      #    println("tolerance reached, ending at iteration ",iter)
	      #    break;
	      # end
    end
    return m, cost
end

function PCG{I<:Int64,C<:Complex128,F<:Float64}(d::Vector{Complex{Float64}},x0::Vector{Complex{F}},
             woz::I,wox::I,wnz::I,wnx::I,irz::Vector{I},irx::Vector{I},M::Mesh,
             W::SparseMatrixCSC{F, I},
             H::Base.SparseArrays.UMFPACK.UmfpackLU{C,I},
             A::Base.SparseArrays.UMFPACK.UmfpackLU{C,I};
             Niter=10, mu=0.)

    J = Float64[]
    x = copy(x0)
    tx= copy(x0)
    s = copy(x0)
    tx= W * x
    t = forwardOperator(tx,irz,irx,woz,wox,wnz,wnx,M,H)
    r = d - t
    cost0 = norm(r)^2 + mu * norm(x)^2
    ta = adjointOperator(r,irz,irx,woz,wox,wnz,wnx,M,A)
    ta = W * ta
    s  = ta - mu * x
    p  = copy(s)
    gamma = norm(s)^2
    for iter = 1 : Niter
        tx = W * p
        q = forwardOperator(tx,irz,irx,woz,wox,wnz,wnx,M,H)
	      delta = norm(q)^2 + mu * norm(p)^2
	      alpha = gamma/delta
	      x = x + alpha * p
	      r = r - alpha * q
        cost = norm(r)^2 + mu * norm(x)^2
        tmp = cost / cost0
        println("iteration $iter, cost $tmp")
	      s = adjointOperator(r,irz,irx,woz,wox,wnz,wnx,M,A)
        s = W * s
	      s = s - mu * x
	      gamma0 = copy(gamma)
        gamma  = norm(s)^2
        beta   = gamma / gamma0
        p      = s + beta * p
        push!(J,tmp)
    end
    x = W * x
    return x, J
end

function goodPass(x::Vector{Complex{Float64}})
    tmp = zeros(Float64, length(x))
    for i = 1 : length(x)
        tmp[i] = sqrt(abs(x[i]))
    end
    return spdiagm(tmp)
end

function stackInversion(dc::Matrix{Float64}, dt::Float64,
                        L::SparseMatrixCSC{Float64,Int64}, m::Matrix{Float64},
                        pml::Matrix{Float64}, somm::Matrix{Complex128},
                        woz::Int64, wox::Int64, wnz::Int64, wnx::Int64,
                        irz::Vector{Int64}, irx::Vector{Int64}, M::Mesh,
                        miniw::Int64, maxiw::Int64)
    nt = size(dc,1)
    fdc= fft(dc,1)
    df = 1./(nt*dt)
    loc = zeros(wnz, wnx, maxiw-miniw+1)
    for iw = miniw : maxiw
        omega = 2*pi*(iw-1)*df
        dobs  = fdc[iw, :]
        H = getHelmhotzMtx(L, omega, m, pml, somm); H = lufact(H);
        A = getAdjointMtx( L, omega, m, pml, somm); A = lufact(A);
        lambda = fpower(woz, wox, wnz, wnx, irz, irx, M, H, A, maxiter=60);
        lambda = ceil(lambda, 10); mu = 2.5*lambda
        (tmp, J) = ffista(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A, mu, lambda, maxiter=100);
        tmp = reshape(tmp, wnz, wnx);
        loc[:,:,iw-miniw+1] = abs(tmp[:,:])
        println("finish $iw")
    end
    return loc
end

function getlambda(miniw::Int64,maxiw::Int64,dt::Float64,nt::Int64,
                   woz::Int64,wox::Int64,wnz::Int64,wnx::Int64,
                   irz::Vector{Int64},irx::Vector{Int64},M::Mesh,
                   L::SparseMatrixCSC{Float64,Int64},m::Matrix{Float64},
                   pml::Matrix{Float64},somm::Matrix{Complex128})
    nw = maxiw - miniw+1
    lambda = zeros(nw)
    dw = 2.0*pi/(nt*dt)
    list = Vector{Tuple{Float64,Int64,Int64,Int64,Int64,
                        Vector{Int64},Vector{Int64},Mesh,
                        SparseMatrixCSC{Float64,Int64},
                        Matrix{Float64},Matrix{Float64},Matrix{Complex128}}}(nw)
    for iw = miniw : maxiw
        omega = (iw-1)*dw
        list[iw-miniw+1] = (omega, woz, wox, wnz, wnx, irz, irx, M, L, m, pml, somm)
    end
    lambda = pmap(getEig, list)
    return lambda
end

function getEig(par::Tuple{Float64,Int64,Int64,Int64,Int64,
                           Vector{Int64},Vector{Int64},Mesh,
                           SparseMatrixCSC{Float64,Int64},
                           Matrix{Float64},Matrix{Float64},Matrix{Complex128}})
    omega = par[1]
    woz = par[2]; wox = par[3]; wnz = par[4]; wnx = par[5];
    irz = par[6]; irx = par[7]; M   = par[8]; L   = par[9];
    m   = par[10]; pml = par[11];somm= par[12];
    H = getHelmhotzMtx(L, omega, m, pml, somm); H = lufact(H);
    A = getAdjointMtx( L, omega, m, pml, somm); A = lufact(A);
    lambda = fpower(woz, wox, wnz, wnx, irz, irx, M, H, A, maxiter=60);
    println("finish $omega")
    return lambda
end

function initRemoteChannel(func::Union{Function,Type},pid::Int64,args...)
    return RemoteChannel(()->initChannel(func,args), pid)
end

function initChannel(func::Union{Function,Type},args::Tuple)
    obj = func(args...)
    chan = Channel{typeof(obj)}(1)
    put!(chan,obj)
    return chan
end

type Location
     M   :: Mesh
     m   :: Matrix{Float64}
     irz :: Vector{Int64}
     irx :: Vector{Int64}
     woz :: Int64
     wox :: Int64
     wnz :: Int64
     wnx :: Int64
     omega :: Float64
     dobs :: Vector{Complex128}
     lambda:: Float64
     H :: Base.SparseArrays.UMFPACK.UmfpackLU{Complex128,Int64}
     A :: Base.SparseArrays.UMFPACK.UmfpackLU{Complex128,Int64}
     function Location(M::Mesh, m::Matrix{Float64}, L::SparseMatrixCSC{Float64,Int64},
                       pml::Matrix{Float64}, somm::Matrix{Complex128},
                       irz::Vector{Int64}, irx::Vector{Int64}, woz::Int64,
                       wox::Int64, wnz::Int64, wnx::Int64, iw::Int64, dobs::Vector{Complex128})
         omega = 2.*pi*(iw-1) / (M.nt * M.dt)
         H = getHelmhotzMtx(L, omega, m, pml, somm); H = lufact(H);
         A = getAdjointMtx( L, omega, m, pml, somm); A = lufact(A);
         lambda = fpower(woz, wox, wnz, wnx, irz, irx, M, H, A)
         lambda = ceil(lambda, 10)
         println("set up $iw")
         return new(M, m, irz, irx, woz, wox, wnz, wnx, omega, dobs, lambda, H, A)
     end
end

function freqHost(miniw::Int64, maxiw::Int64)
    np = nworkers()
    fid = Vector{Vector{Int64}}(np)
    nw = maxiw - miniw + 1
    nj = floor(Int64, nw/np)
    res= rem(nw, np)
    for i = 1 : np
        il = miniw + (i-1)*nj
        iu = il + nj - 1
        fid[i] = collect(il:iu)
    end
    for i = 1 : res
        push!(fid[i], i+nj*np+miniw-1)
    end
    return fid
end

function setupLocation(miniw::Int64, maxiw::Int64, M::Mesh, m::Matrix{Float64},
                       L::SparseMatrixCSC{Float64,Int64}, pml::Matrix{Float64},
                       somm::Matrix{Complex128}, irz::Vector{Int64}, irx::Vector{Int64},
                       woz::Int64, wox::Int64, wnz::Int64, wnx::Int64, dc::Matrix{Float64})
    fid = freqHost(miniw, maxiw)
    nw  = maxiw - miniw + 1
    fdc = fft(dc, 1)
    nq  = length(fid)
    workerList = workers()
    rr = Vector{RemoteChannel}(nw)
    @sync begin
          for p = workerList
              @async begin
                     for i = 1 : length(fid[p-1])
                         iw = fid[p-1][i]
                         dobs = fdc[iw,:]
                         rr[fid[p-1][i]-miniw+1] = initRemoteChannel(Location, p, M, m,
                                                   L, pml, somm, irz, irx, woz, wox, wnz,
                                                   wnx, iw, dobs)
                     end
              end
          end
    end
    return rr
end

function updateDobs(rr::Vector{RemoteChannel}, dobs::Matrix{Float64}, miniw::Int64, maxiw::Int64)
    fid = freqHost(miniw, maxiw)
    fdc = fft(dobs, 1)
    workerList = workers()
    @sync begin
          for p = workerList
              @async begin
                     for i = 1 : length(fid[p-1])
                         iw  = fid[p-1][i]
                         idx = fid[p-1][i]-miniw+1
                         slice = fdc[iw,:]
                         remotecall_fetch(updateDobs, p, rr[idx], slice)
                     end
              end
          end
    end
    return nothing
end

function updateDobs(rc::RemoteChannel, slice::Vector{Complex128})
    loc = take!(rc)
    loc.dobs = slice
    put!(rc, loc)
    return nothing
end

function remoteFista(rr::RemoteChannel{Channel{AcousticWave.Location}}, mu, maxiter)
    loc = fetch(rr)
    mu1 = loc.lambda * mu
    freq = loc.omega / (2*pi)
    (x,J) = ffista(loc.dobs, loc.woz, loc.wox, loc.wnz, loc.wnx,
                   loc.irz, loc.irx, loc.M, loc.H, loc.A, mu1, loc.lambda)
    println("$freq")
    return abs(x)
end

function parLocation(rr::Vector{RemoteChannel}, wnz::Int64, wnx::Int64, miniw::Int64, maxiw::Int64; mu=2.5, maxiter=100)
    nw = length(rr)
    workerList = workers()
    amp = zeros(wnz*wnx, nw)
    fid = freqHost(miniw, maxiw)
    @sync begin
          for p = workerList
              @async begin
                     for i = 1 : length(fid[p-1])
                         idx = fid[p-1][i] - miniw + 1
                         amp[:,idx] = remotecall_fetch(remoteFista, p, rr[idx], mu, maxiter)
                     end
              end
          end
    end
    return amp
end

# function test(rr::Vector{RemoteChannel})
#     n = length(rr)
#     lambda = Vector{Float64}(n)
#     @sync begin
#           for i = 1 : n
#               @async begin
#                      ip = rr[i].where
#                      lambda[i] = remotecall_fetch(()->fetch(rr[i]).lambda, ip)
#               end
#           end
#     end
#     return lambda
# end
#
# function test1(rr::Vector{RemoteChannel})
#     n = length(rr)
#     lambda = Vector{Location}(n)
#     @sync begin
#           for i = 1 : n
#               @async begin
#                      ip = rr[i].where
#                      lambda[i] = remotecall_fetch(()->fetch(rr[i]).lambda, ip)
#               end
#           end
#     end
#     return lambda
# end

# function IRLS{I<:Int64, F<:Float64}(dobs::Vector{Complex{Float64}},
#              woz::I,wox::I,wnz::I,wnx::I,irz::Vector{I},irx::Vector{I},M::Mesh,
#              H::Base.SparseArrays.UMFPACK.UmfpackLU{Complex{F},I},
#              A::Base.SparseArrays.UMFPACK.UmfpackLU{Complex{F},I};
#              maxiter=5, Niter=10, mu=0.000001)
#      (m, J) = CG(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A, maxiter=maxiter)
#      for i = 1 : Niter
#          P = goodPass(m)
#          (m, J) = PCG(dobs,woz,wox,wnz,wnx,irz,irx,M,P,H,A,maxiter=maxiter)
#      end
#      return reshape(m, wnz, wnx)
# end

# dot product test of source location
# m = rand(Complex128, wnz*wnx);
# d = forwardOperator(m, woz, wox, wnz, wnx, irz, irx, M, H);
#
# d1 = rand(Complex128, length(irz));
# m1 = adjointOperator(d1, woz, wox, wnz, wnx, irz, irx, M, A);
# tmp = dot(m, m1)
# tmp1= dot(d, d1)
