using SeisPlot, SeisAcoustic, LinearAlgebra

# working directory
work_dir = joinpath(homedir(), "Desktop/tmp_LSRTM");

# size of model
nz = 101; nx = 301;

# 2-layers velocity model
vel = 3000 * ones(nz, nx);  # m/s
# vel[51:end,:] .= 3500;

# constant density model
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = false;

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;

# precision
data_format=Float64;
order=3

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# initialize a source
# isz = 2; isx = 150;
# src = Source(isz, isx, params; amp=100000, fdom=15);

isx = collect(5:10:295); ns=length(isx);
isz = 2   * ones(ns);
ot  = 0.5 * rand(ns);
pol = rand(Bool, ns);
amp = 100000 * ones(ns);
for i = 1 : ns
    if !pol[i]
       amp[i] = - amp[i]
    end
end
srcs = get_multi_sources(isz, isx, params; amp=amp, ot=ot, fdom=5);

# generate observed data
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
dobs= Recordings(irz, irx, params);

# generate recordings and boundary value
path_bnd = joinpath(work_dir, "bnd.bin")
path_wfd = joinpath(work_dir, "wfd.bin")
path_sws = joinpath(work_dir, "sws.bin")

@time multi_step_forward!(dobs, srcs, params; path_bnd=path_bnd, path_wfd=path_wfd, path_sws=path_sws);
@time pre = get_sourceside_wavefield(srcs, params);

(hdr, s1) = read_RSdata(path_sws);
s2 = zeros(nz, nx)
for i2 = 1 : nx
    for i1 = 1 : nz
        s2[i1,i2] = dot(pre[i1,i2,:], pre[i1,i2,:])
    end
end
norm(s1-s2) / norm(s1)

SeisPlotTX(dobs.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

a = b
# ==============================================================================
#                      test sources-side wavefield
# ==============================================================================
using SeisPlot, SeisAcoustic

#             read physical model
work_dir = joinpath(homedir(), "Desktop/tmp_LSRTM");

path_rho = joinpath(work_dir, "physical_model/rho.rsf");
path_vel = joinpath(work_dir, "physical_model/vel.rsf");

# read the physical model
(hdr_rho, rho) = read_RSdata(path_rho);
(hdr_vel, vel) = read_RSdata(path_vel);

# # cropped model for imaging
z_range = 1:2:300;
x_range = 4000:2:4600;
vel = vel[z_range, x_range];
rho = rho[z_range, x_range];
vel = model_smooth(vel, 10);

# SeisPlotTX(vel, hbox=1.5*2, wbox=3*2, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
# SeisPlotTX(rho, hbox=1.5*2, wbox=3*2, cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0008; tmax = 1.5;

# top boundary condition
free_surface = false;
data_format  = Float32;
order        = 5;

# get the source-side wavefield boundary
rho1 = copy(rho); rho1 .= minimum(rho);  # homogenous density model
params = TdParams(rho1, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# a single source
# isz = 2; isx = 50;
# src = Source(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# vector of source
isx = collect(30:70:params.nx-30); isz = 2*ones(length(isx));
src = get_multi_sources(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

dir_sourceside = joinpath(work_dir, "sourceside");
get_wavefield_bound(dir_sourceside, src, params; remove_flag=false);

(hdr1, d1) = read_RSdata(joinpath(dir_sourceside, "normalization.rsf"));
SeisPlotTX(d1, dy=dt, cmap="rainbow", wbox=6, hbox=3);

(h1, s1) = read_RSdata(joinpath(dir_sourceside, "strength/strength_1.rsf"));
(h2, s2) = read_RSdata(joinpath(dir_sourceside, "strength/strength_2.rsf"));
(h3, s3) = read_RSdata(joinpath(dir_sourceside, "strength/strength_3.rsf"));
(h4, s4) = read_RSdata(joinpath(dir_sourceside, "strength/strength_4.rsf"));
SeisPlotTX(s1, dy=dt, cmap="rainbow", wbox=6, hbox=3);
SeisPlotTX(s1+s2, dy=dt, cmap="rainbow", wbox=6, hbox=3);
SeisPlotTX(s1+s2+s3, dy=dt, cmap="rainbow", wbox=6, hbox=3);
SeisPlotTX(s1+s2+s3+s4, dy=dt, cmap="rainbow", wbox=6, hbox=3);

# apply the adjoint operator
path_m = joinpath(work_dir, "model_space/m.rsf");
dir_rec= joinpath(work_dir, "data_space/observation");
dir_sourceside = joinpath(work_dir, "sourceside");
get_born_adjoint(path_m, dir_rec, dir_sourceside, params; path_)
