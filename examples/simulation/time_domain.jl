using SeisPlot, SeisAcoustic, LinearAlgebra, DSP

# homogeneous velocity and density model
nz = 100; nx = 201;
vel = 2800 * ones(nz, nx);  # m/s
vel[31:end,:] .= 3100;  # m/s
vel[61:end,:] .= 3500;  # m/s
rho = 2.0 * ones(nz, nx);  # kg/m^3

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

# dot-product test (linear operator
irx = collect(1:2:params.nx); irz = 2 * ones(length(irx));
rec1= multi_step_forward!(irz, irx, src, params);

# band limited random recordings
w    = ricker(10.0, params.dt);
hl   = floor(Int64, length(w)/2);
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

# ==============================================================================
#                          simultaneouse source
# ==============================================================================
isx = collect(5:40:params.nx); ns=length(isx); isz = 2*ones(ns);
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

irx = collect(1:2:params.nx); irz = 2 * ones(length(irx));
rec1= multi_step_forward!(irz, irx, src, params);

# band limited random recordings
w    = ricker(10.0, params.dt);
hl   = floor(Int64, length(w)/2);
rec2 = Recordings(irz, irx, params);
idx_o= [1]
for i = 1 : rec2.nr
    tmp = conv(randn(params.nt)*1000, w)
    copyto!(rec2.p, idx_o[1], tmp, hl+1, params.nt)
    idx_o[1] = idx_o[1] + params.nt
end
p1 = multi_step_adjoint!(rec2, src, params);

tmp1 =[0.0];
for i = 1 : ns
    tmp1[1] = tmp1[1] + dot(p1[i], src[i].p)
end
tmp2 = dot(vec(rec2.p), vec(rec1.p))
(tmp1[1]-tmp2) / tmp2

# save the adjoint wavefield
path_spt = joinpath(homedir(), "Desktop/snapshot.rsb" );
path_wfd = joinpath(homedir(), "Desktop/wavefield.rsb");
path_pre = joinpath(homedir(), "Desktop/pressure.rsb" );

multi_step_adjoint_test(rec, params; path_spt=path_spt, path_wfd=path_wfd,
                    path_pre = path_pre, interval=1);

(hdr_s, s) = read_RSdata(path_spt);
(hdr_w, w) = read_RSdata(path_wfd);
(hdr_p, p) = read_RSdata(path_pre);
i = 20;
SeisPlotTX(p[:,:,i], pclip=98, cmap="seismic", wbox=10, hbox=5);

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
