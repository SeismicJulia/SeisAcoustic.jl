using SeisPlot, SeisAcoustic, LinearAlgebra, DSP

# homogeneous velocity and density model
nz = 101; nx = 301;
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# number of PML layers
npml = 20;

# top boundary condition
free_surface = true;   #(pml or free_surface)

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 3.0;  # use second as unit

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);

# shows the default value for the keyword parameters
# data_format = (Float32 or Float64)
# fd_flag     = ("taylor" or "ls")
# order       = 2 - 10 if we use "ls" to compute the FD coefficients
#             = 2 - n  if we use "taylor" expansion to compute FD coefficients

# ==============================================================================
#                          single source simulation
# ==============================================================================
src = Source(2, 150, params; ot=0.0, fdom=15.0,
      type_flag="ricker", amp=100000, location_flag="index");

# initialize recordings
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));

# forward modeling of simultaneous sources
rec       = multi_step_forward!(irz, irx, src, params);
path_shot = multi_step_forward!(irz, irx, src, params; path_shot="/Users/wenlei/Desktop/test.bin");
rec1      = read_recordings(path_shot);
recordings_isequal(rec, rec1)

# save the boundary value and last wavefield
path_spt = joinpath(homedir(), "Desktop/snapshot.rsb" );
path_wfd = joinpath(homedir(), "Desktop/wavefield.rsb");
path_pre = joinpath(homedir(), "Desktop/pressure.rsb" );
path_bnd = joinpath(homedir(), "Desktop/boundary.rsb" );
path_lwfd= joinpath(homedir(), "Desktop/last_wfd.rsb" );
path_sws = joinpath(homedir(), "Desktop/strength.rsb" );

multi_step_forward!(src, params; path_spt=path_spt, path_wfd=path_wfd, path_pre=path_pre, interval=1,
                    path_bnd=path_bnd, path_lwfd=path_lwfd, path_sws=path_sws);

# save the pressure field
(hdr, p0) = read_RSdata(path_pre);

# reconstruct the pressure field forward
p1 = pressure_reconstruct_forward(path_bnd, src, params);

# reconstruct the pressure field backward
p2 = pressure_reconstruct_backward(path_bnd, path_lwfd, src, params);

# testing wavefield reconstruction forward or backward
SeisPlotTX(p0[:,:,500], pclip=98, cmap="seismic", wbox=6, hbox=2);
SeisPlotTX(p1[:,:,500], pclip=98, cmap="seismic", wbox=6, hbox=2);
SeisPlotTX(p2[:,:,500], pclip=98, cmap="seismic", wbox=6, hbox=2);
norm(p0-p1) / norm(p0)
norm(p0-p2) / norm(p0)

# ==============================================================================
#                          simultaneouse source
# ==============================================================================
isx = collect(5:60:295); ns=length(isx); isz = 2*ones(ns);
ot  = 0.5*rand(ns);
src = get_multi_sources(isz, isx, params; amp=100000, ot=ot, fdom=15);

# forward modeling of simultaneous sources
rec       = multi_step_forward!(irz, irx, src, params);
path_shot = multi_step_forward!(irz, irx, src, params; path_shot="/Users/wenlei/Desktop/test.bin");
rec1      = read_recordings(path_shot);
recordings_isequal(rec, rec1)

# save the boundary value and last wavefield
path_spt = joinpath(homedir(), "Desktop/snapshot.rsb" );
path_wfd = joinpath(homedir(), "Desktop/wavefield.rsb");
path_pre = joinpath(homedir(), "Desktop/pressure.rsb" );
path_bnd = joinpath(homedir(), "Desktop/boundary.rsb" );
path_lwfd= joinpath(homedir(), "Desktop/last_wfd.rsb" );
path_sws = joinpath(homedir(), "Desktop/strength.rsb" );

multi_step_forward!(src, params; path_spt=path_spt, path_wfd=path_wfd, path_pre=path_pre, interval=1,
                    path_bnd=path_bnd, path_lwfd=path_lwfd, path_sws=path_sws);

# save the pressure field
(hdr, p0) = read_RSdata(path_pre);

# reconstruct the pressure field forward
p1 = pressure_reconstruct_forward(path_bnd, src, params);

# reconstruct the pressure field backward
p2 = pressure_reconstruct_backward(path_bnd, path_lwfd, src, params);

# testing wavefield reconstruction forward or backward
SeisPlotTX(p0[:,:,500], pclip=98, cmap="seismic", wbox=6, hbox=2);
SeisPlotTX(p1[:,:,500], pclip=98, cmap="seismic", wbox=6, hbox=2);
SeisPlotTX(p2[:,:,500], pclip=98, cmap="seismic", wbox=6, hbox=2);
norm(p0-p1) / norm(p0)
norm(p0-p2) / norm(p0)

# ==============================================================================
#                   dot-product test forward and adjoint operator
# ==============================================================================
# test the one-step adjoint operator
spt1_f = Snapshot(params);
spt2_f = Snapshot(params);
spt1_b = Snapshot(params);
spt2_b = Snapshot(params);

# initialize spt1_f with random number
for ix = 1 : params.Nx
    amp = 1.0
    col_idx = (ix-1) * params.Nz
    for iz = 1 : params.Nz
        idx= col_idx + iz
        spt1_f.vz[idx] = amp * randn(); spt2_b.vz[idx] = amp * randn()
        spt1_f.vx[idx] = amp * randn(); spt2_b.vx[idx] = amp * randn()
        spt1_f.pz[idx] = amp * randn(); spt2_b.pz[idx] = amp * randn()
        spt1_f.px[idx] = amp * randn(); spt2_b.px[idx] = amp * randn()
    end
end

# temporary variables
tmp    = zeros(params.Nz * params.Nx);
tmp_z1 = zeros(params.data_format, params.Nz);
tmp_z2 = zeros(params.data_format, params.Nz);
tmp_x1 = zeros(params.data_format, params.Nx);
tmp_x2 = zeros(params.data_format, params.Nx);

# nt-step forward
nt = 200
for it = 1 : nt
    one_step_forward!(spt2_f, spt1_f, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2);
    copy_snapshot!(spt1_f, spt2_f);
end

# nt-step adjoint
for it = 1 : nt
    one_step_adjoint!(spt1_b, spt2_b, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2);
    copy_snapshot!(spt2_b, spt1_b);
end

# inner product
tmp1 = (dot(spt1_f.vz, spt1_b.vz) + dot(spt1_f.vx, spt1_b.vx)
      + dot(spt1_f.pz, spt1_b.pz) + dot(spt1_f.px, spt1_b.px))

tmp2 = (dot(spt2_f.vz, spt2_b.vz) + dot(spt2_f.vx, spt2_b.vx)
      + dot(spt2_f.pz, spt2_b.pz) + dot(spt2_f.px, spt2_b.px))

(tmp1-tmp2) / tmp1


# ==============================================================================
#                   dot-product test forward and adjoint operator
# ==============================================================================
irx = collect(1:2:params.nx); irz = 2 * ones(length(irx));
rec1= multi_step_forward!(irz, irx, src, params);

# band limited random recordings
w    = ricker(10.0, params.dt)
hl   = floor(Int64, length(w)/2)
rec2 = Recordings(irz, irx, params);
idx_o= [1]
for i = 1 : rec2.nr
    tmp = conv(randn(params.nt)*1000, w)
    copyto!(rec2.p, idx_o[1], tmp, hl+1, params.nt)
    idx_o[1] = idx_o[1] + params.nt
end
p1 = multi_step_adjoint!(rec2, src, params);

tmp1 = dot(p1, src.p)
tmp2 = dot(vec(rec2.p), vec(rec1.p))
(tmp1-tmp2) / tmp1

# using SeisPlot, SeisAcoustic
#
# # homogeneous velocity and density model
# path_vel = joinpath(homedir(), "Desktop/vp.bin")
# fid = open(path_vel)
# vel = zeros(Float32, 750*150); read!(fid, vel); close(fid);
# vel = reshape(vel, 750, 150); vel = permutedims(vel, [2,1]);
#
# path_rho = joinpath(homedir(), "Desktop/rho.bin");
# fid = open(path_rho)
# rho = zeros(Float32, 750*150); read!(fid, rho); close(fid);
# rho = reshape(rho, 750, 150); rho = permutedims(rho, [2,1]);
#
# # number of PML layers
# npml = 20;
#
# # top boundary condition
# free_surface = true;
#
# # vertical and horizontal grid size
# dz = 10; dx = 10;
#
# # time step size and maximum modelling length
# dt = 0.0015; tmax = 2.25;  # use second as unit
#
# # organize these parameters into a structure
# params2 = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
#          data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);
# params6 = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
#          data_format=Float64, fd_flag="taylor", order=6, npml=20, apml=900.);
#
# # initialize a source
# src = Source(2, 375, params; ot=0.0, fdom=10.0,
#       type_flag="ricker", amp=100000, location_flag="index");
#
#
# # initialize recordings
# irx = collect(1:2:params.nx);
# irz = 2 * ones(length(irx));
# rec2 = Recordings(irz, irx, params);
# rec6 = Recordings(irz, irx, params);
#
#
# # forward modeling of simultaneous sources
# multi_step_forward!(rec2, src , params2);
# multi_step_forward!(rec6, src , params6);
# SeisPlotTX(rec2.p, wbox=10, hbox=10, cmap="seismic");
# SeisPlotTX(rec6.p, wbox=10, hbox=10, cmap="seismic");
# SeisPlotTX(rec6.p-rec2.p, wbox=10, hbox=10, cmap="seismic");
#
# tmp = hcat(rec2.p, rec6.p, rec6.p-rec2.p);
# SeisPlotTX(tmp, wbox=30, hbox=10, cmap="seismic");
#
# # save the pressure field
# path_pre2 = joinpath(homedir(), "Desktop/pressure2.rsb");
# path_pre6 = joinpath(homedir(), "Desktop/pressure6.rsb");
#
# # multi_step_forward!(path_pre, src, params; save_flag="pressure");
# multi_step_forward!(path_pre2, src, params2; save_flag="pressure");
# multi_step_forward!(path_pre6, src, params6; save_flag="pressure");
# (hdr, p0) = read_RSdata(path_pre6);
#
# it = 300; SeisPlotTX(p0[:,:,it], wbox=15, hbox=3, cmap="seismic", name="/Users/wenlei/Desktop/spt300.pdf");
# it = 600; SeisPlotTX(p0[:,:,it], wbox=15, hbox=3, cmap="seismic", name="/Users/wenlei/Desktop/spt600.pdf");
# it = 900; SeisPlotTX(p0[:,:,it], wbox=15, hbox=3, cmap="seismic", name="/Users/wenlei/Desktop/spt900.pdf");
