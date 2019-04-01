using AcousticWave

# generate synthetic data in time domain
function wraperFdModeling()
    nz = 150; nx = 300; ext = 20; iflag = 1;
    dz = 8.; dx = 8.; dt = 1.e-3; f0=20.0; tmax = 4.5;
    v  = 3000. * ones(nz, nx);

    # first stage
    isz0 = collect(55:2:95)            ; tmpz = collect(55:2:95)            ; ot0 = vcat(collect(0.5:-0.05:0.), collect(0.05:0.05:0.5));
    isx0 = 197*ones(Int64,length(isz0)); tmpx = 202*ones(Int64,length(tmpz)); tt  = vcat(collect(0.5:-0.05:0.), collect(0.05:0.05:0.5));
    isz0 = vcat(isz0, tmpz); isx0 = vcat(isx0, tmpx); ot0 = vcat(ot0, tt);

    # second stage
    isz1 = collect(65:2:85)            ; tmpz = collect(65:2:85)            ; ot1 = vcat(collect(0.25:-0.05:0.), collect(0.05:0.05:0.25));
    isx1 = 148*ones(Int64,length(isz1)); tmpx = 153*ones(Int64,length(tmpz)); tt  = vcat(collect(0.25:-0.05:0.), collect(0.05:0.05:0.25));
    isz1 = vcat(isz1, tmpz); isx1 = vcat(isx1, tmpx); ot1 = vcat(ot1, tt);
    ot1  = ot1 + 1.5;

    # third stage
    isz2 = collect(65:2:85)           ; tmpz = collect(55:2:95)            ; ot2 = vcat(collect(0.25:-0.05:0.), collect(0.05:0.05:0.5));
    isx2 = 94*ones(Int64,length(isz2)); tmpx = 104*ones(Int64,length(tmpz)); tt  = vcat(collect(0.25:-0.05:0.), collect(0.05:0.05:0.5));
    isz2 = vcat(isz2, tmpz); isx2 = vcat(isx2, tmpx); ot2 = vcat(ot2, tt);
    ot2  = ot2 + 3.0;

    # all the sources
    isz = vcat(isz0, isz1, isz2);
    isx = vcat(isx0, isx1, isx2);
    ot  = vcat(ot0, ot1, ot2);

    amp = 1.0*ones(length(isz));
    srcs= InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
    # discretize spatial derivative operator

    vmax = maximum(v); vmin = minimum(v);
    fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);

    irz = collect(1:2:nz); irx = 5*ones(Int64, length(irz));
    tmpx = collect(1:2:nx); tmpz = 5*ones(Int64, length(tmpx));
    irz = vcat(irz, tmpz); irx = vcat(irx, tmpx);
    # irx = collect(1:1:nx); irz = 5*ones(Int64, length(irx));
    # shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax);
    # return shot.p
    path = "/Users/wenlei/Desktop/wfd.bin"
    MultiStepForward(path, srcs, fidMtx, tmax = tmax)
    return nothing
end
@time d = wraperFdModeling();

path = "/Users/wenlei/Desktop/moive"
path1= "/Users/wenlei/Desktop/wfd.bin"
waveAnim (path, path1)

v = quantile(vec(d),0.95)
imshow(d, cmap="gray", aspect="auto", vmax=v, vmin=-v)

ot = 0.1; nt=1200; dt=0.001; df = 1./(nt*dt);
dc0 = intercept(d, ot, dt, nt);
imshow(dc0, cmap="gray", aspect="auto", vmax=v, vmin=-v);

ot = 1.55; dt=0.001; df = 1./(nt*dt);
dc1= intercept(d, ot, dt, nt);
figure(); imshow(dc1, cmap="gray", aspect="auto", vmax=v, vmin=-v);

ot = 3.; dt=0.001; df = 1./(nt*dt);
dc2= intercept(d, ot, dt, nt);
figure(); imshow(dc2, cmap="gray", aspect="auto", vmax=v, vmin=-v);
# # # Source location by stacking FISTA result
# # # ========================================
nz=150; nx=300; ext=50; iflag=1;
dz=0.008; dx=0.008;
M  = Mesh(nt, nz, nx, ext, iflag, dt, dz, dx);
# model parameters
v = 3. * ones(nz, nx); m = 1.0 ./ (v.^2); m = modelExpand(m, M);
# boundary conditions
amp = -1.5*maximum(v)./(ext*dz)*log(0.001); pml= getABC(M, amp);
somm = getSommerfeldBC(M, m);
# laplace operator
L = getLaplacian(M);
# receivers' location
irz  = collect(1:2:nz); irx = 5*ones(Int64, length(irz));
tmpx = collect(1:2:nx); tmpz= 5*ones(Int64, length(tmpx));
irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
# irx = collect(1:1:nx); irz = 5*ones(Int64, length(irx));
# window size and location
woz=10; wox=10; wnz=140; wnx=290; nt=1024;
miniw = 5; maxiw = 55;

# set up location parameters
rr = setupLocation(miniw, maxiw, M, m, L, pml, somm,
                   irz, irx, woz, wox, wnz, wnx, dc0)

# begin source location
# first stage
updateDobs(rr, dc0, miniw, maxiw);
amp0 = parLocation(rr, wnz, wnx, miniw, maxiw; mu=3.0);
loc0 = reshape(sum(amp0,2), wnz, wnx);

# second stage
updateDobs(rr, dc1, miniw, maxiw);
amp1 = parLocation(rr, wnz, wnx, miniw, maxiw; mu=3.0)
loc1 = reshape(sum(amp1,2), wnz, wnx);

# third stage
updateDobs(rr, dc2, miniw, maxiw);
amp2 = parLocation(rr, wnz, wnx, miniw, maxiw; mu=3.0)
loc2 = reshape(sum(amp2,2), wnz, wnx);

Nz=nz+2*ext; Nx=nx+2*ext;
v = quantile(vec(loc0), 1.0)
s0 = zeros(Nz, Nx); s0[ext+woz+1:ext+woz+wnz, ext+wox+1:ext+wox+wnx]=loc0[:,:]
imshow(s0, vmax=v, vmin=0);

v = quantile(vec(loc1), 1.0)
s1 = zeros(Nz, Nx); s1[ext+woz+1:ext+woz+wnz, ext+wox+1:ext+wox+wnx]=loc1[:,:]
figure(); imshow(s1, vmax=v, vmin=0);

v = quantile(vec(loc2), 1.0)
s2 = zeros(Nz, Nx); s2[ext+woz+1:ext+woz+wnz, ext+wox+1:ext+wox+wnx]=loc2[:,:]
figure(); imshow(s2, vmax=v, vmin=0);

function hybrid(amp::Matrix{Float64}, n::Int64)
    (nr, nc) = size(amp)
    ns = floor(Int64, nc/n)
    if ns * n < nc
       n = n + 1
    end
    tmp = zeros(nr, n)
    for i = 1 : n
        icl = (i-1)*ns+1
        icu = icl + ns
        if icu > nc
           icu = nc
        end
        for j = icl : icu
            tmp[:,i] = tmp[:,i] + amp[:,j]
        end
    end
    return tmp
end




# test the result of observaed data updating
# fdc1 = fft(dc1, 1); idx = 47;
# iw = idx + miniw - 1
# spl = fetch(rr[idx])
# norm(spl.dobs-fdc1[iw,:])

# loc = Location(M, m, L, pml, somm, irz, irx, woz, wox, wnz, wnx, iw, dc[20,:])
# lambda = getlambda(miniw,maxiw,dt,nt,woz,wox,wnz,wnx,irz,irx,M,L,m,pml,somm)

# iw = 30; df=0.997; omega = 2*pi*(iw-1)*df
# # fdc = fft(dc,1); dobs = fdc[iw, :];
# H = getHelmhotzMtx(L, omega, m, pml, somm); H = lufact(H);
# A = getAdjointMtx( L, omega, m, pml, somm); A = lufact(A);
#
# # estimate the maximum eigenvalue
# lambda = fpower(woz, wox, wnz, wnx, irz, irx, M, H, A, maxiter=60);
#
# lambda = 4.096e-6; mu = 2e-5;
# (x40, J) = ffista(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A, mu, lambda, maxiter=100);
# tmp40 = reshape(abs(x40), wnz, wnx); imshow(tmp40);
#
# miniw = 10; maxiw= 55;
# loc = stackInversion(dc, dt, L, m, pml, somm,
#                      woz, wox, wnz, wnx, irz, irx, M, miniw, maxiw)
# path = "/Users/wenyue/Desktop/loc1.bin"
# (n1,n2,n3) = size(loc);
# fid = open(path,"w")
# write(fid,n1,n2,n3);
# write(fid,vec(loc));
# close(fid);
#
# tmp = squeeze(sum(abs(loc),3),3);
# v = quantile(vec(tmp), 0.99)
# imshow(tmp, vmax=v, vmin=0)
#
#
#
# x = collect(1:6);
# y = collect(1:6);
# rr = Vector{RemoteChannel}(6)
# assignment = Vector{Vector{Int64}}(nworkers())
# assignment[1] = collect(1:2)
# assignment[2] = collect(3:4)
# assignment[3] = collect(5:6)
# @sync begin
#       for (idx, ip) in enumerate(workers())
#           @async begin
#                  for i = 1 : length(assignment[idx])
#                      subidx = assignment[idx][i]
#                      rr[subidx] = initRemoteChannel(point, ip, x[subidx], y[subidx])
#                  end
#           end
#       end
# end

# # Source location by simply stacking adjoint wavefield
# # ====================================================
# nz=150; nx=300; ext=50; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=2048;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# amp=500.; pml= getABC(M, amp);
# m  = 1/9 * ones(nz, nx); m = modelExpand(m, M);
# L  = getLaplacian(M);
# irz  = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:1:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
#
# #shift the recordings
# dc = intercept(d, 0.4, 0.001, 1024)
# SeisPlot(dc, pclip=95, cmap="gray")
#
# miniw=10; maxiw=60; woz=10; wox=10; wnz=140; wnx=290;
# loc = adjoint_stack(dc, dt, miniw, maxiw, L, m, pml, somm, woz, wox, wnz, wnx, irz, irx, M)

# # # Source location by CG method
# # # ====================================================
# nz=150; nx=300; ext=30; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=2048;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# amp=500.; pml= getABC(M, amp);
# m  = 1/9 * ones(nz, nx); m = modelExpand(m, M);
# L  = getLaplacian(M);
#
# irz  = collect(1:10:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:20:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
#
# #shift the recordings
# dc = intercept(d, 0.300, 0.001, 1200)
# SeisPlot(dc, pclip=95, cmap="gray")
#
# woz=10; wox=10; wnz=140; wnx=290; maxiw=110;
# loc = CG_stack(dc, dt, maxiw, L, m, pml, M,
#                woz, wox, wnz, wnx, irz, irx)
#
# # Source location by IRLS method
# # ====================================================
# nz=150; nx=300; ext=50; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=2048;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# amp=500.; pml= getABC(M, amp);
# m  = 1/9 * ones(nz, nx); m = modelExpand(m, M);
# L  = getLaplacian(M);
#
# irz  = collect(1:10:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:20:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
#
# #shift the recordings
# dc = intercept(d, 0.300, 0.001, 1200)
# SeisPlot(dc, pclip=95, cmap="gray")
#
# woz=10; wox=10; wnz=140; wnx=290; iw=30;
# df = 1/(dt*size(dc,1)); f = (iw-1)*df; omega= 2*pi*f;
# H = getHelmhotzMtx(L, omega, m, pml); H = lufact(H);
# A = getAdjointMtx( L, omega, m, pml); A = lufact(A);
#
# # dc = fft(dc,1); dobs = dc[iw,:];
# (loc2, J) = CG(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A, Niter=10, mu=0.)
# for i = 1 : 5
# P = goodPass(loc2)
# (loc2, J) = PCG(dobs,loc1,woz,wox,wnz,wnx,irz,irx,M,P,H,A,mu=1e-16)
# end
# imshow(reshape(abs(loc2), wnz, wnx))

# modeling in frequency domai
# =================================================
# using PyPlot, AcousticWave
# nz=150; nx=300; ext=50  ; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=1024;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# L  = getLaplacian(M);
# # boundary coefficients
# amp = -1.5*3./(ext*dz)*log(0.001); pml= getABC(M, amp);
# # model parameters
# v = 3. * ones(nz, nx);  m = 1./(v.^2);
# m = modelExpand(m, M);
# somm = getSommerfeldBC(M, m);
# # sources
# isz = [50]; isx = [100]; f0 = 20.; t0=[0.];
# src = FSource(isz, isx, f0, t0, dt, nf);
# # receives
# irz  = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:1:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
# # modeling
# maxiw = 85;
# wfd = getDataTime(src, M, pml, somm, L, maxiw)
# # plot the result
# wfd = real(wfd);
# idx = 330; v = maximum(abs(wfd[idx,:,:]));
# figure(); imshow(wfd[idx,:,:], cmap="seismic", vmax=v, vmin=-v)
# savefig("/Users/wenyue/Desktop/wfd.pdf")
# figure(); plot(wfd[idx,:,150])
# savefig("/Users/wenyue/Desktop/trace.pdf")

# function wraperFdModeling()
#     nz = 100; nx = 200; ext = 20; iflag = 1;
#     dz = 8.; dx = 8.; dt = 1.e-3; f0=20.0; tmax = 1.047;
#     v  = 3000. * ones(nz, nx);
#     # initialize source term
#     isz = [50];
#     isx = [100];
#      ot = [0.0];
#     amp = [1.0];
#     srcs= InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
#     # discretize spatial derivative operator
#     vmax = maximum(v); vmin = minimum(v);
#     fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);
#     irz = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
#     tmpx = collect(1:1:nx); tmpz = 5*ones(Int64, length(tmpx));
#     irz = vcat(irz, tmpz); irx = vcat(irx, tmpx);
#     path = "/Users/wenyue/Desktop/wfd.bin"
#     MultiStepForward(path, srcs, fidMtx, tmax=tmax, wtype="p");
# end
# wraperFdModeling()
# spt = readStress(path, idx);
# v = maximum(abs(spt.p));
# imshow(spt.p, cmap="seismic", vmax=v, vmin=-v)
