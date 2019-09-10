using SeisPlot, SeisAcoustic

# ==============================================================================
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

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# homogenous model to simulate direct arrival
rho1 = copy(rho); rho1 .= minimum(rho);
vel1 = copy(vel); vel1 .= vel[1];

params_homo = TdParams(rho1, vel1, free_surface, dz, dx, dt, tmax;
                       data_format=data_format, order=order);

# a single source
# isz = 2; isx = 50;
# src = Source(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# vector of source
isx = collect(30:70:params.nx-30); isz = 2*ones(length(isx));
src = get_multi_sources(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# receiver location
irx = collect(1:2:params.nx);
irz = 2 * ones(Int64, length(irx));

# do simulation
dir_obs = joinpath(work_dir, "data_space/observation")
get_reflections(dir_obs, irz, irx, src, params, params_homo)

# plotting the data
dobs = read_recordings(joinpath(dir_obs, "reflection_4.bin"));
SeisPlotTX(dobs.p, dy=dt, cmap="gray", wbox=6, hbox=6);


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

(h1, s1) = read_RSdata(joinpath(path_sourceside, "strength_1.rsf"));
(h2, s2) = read_RSdata(joinpath(path_sourceside, "strength_2.rsf"));
(h3, s3) = read_RSdata(joinpath(path_sourceside, "strength_3.rsf"));
(h4, s4) = read_RSdata(joinpath(path_sourceside, "strength_4.rsf"));
SeisPlotTX(s1, dy=dt, cmap="rainbow", wbox=6, hbox=3);
SeisPlotTX(s1+s2, dy=dt, cmap="rainbow", wbox=6, hbox=3);
SeisPlotTX(s1+s2+s3, dy=dt, cmap="rainbow", wbox=6, hbox=3);
SeisPlotTX(s1+s2+s3+s4, dy=dt, cmap="rainbow", wbox=6, hbox=3);

# apply the adjoint operator
path_m = joinpath(work_dir, "model_space/m.rsf")
dir_rec= joinpath(work_dir, "data_space/observation")
dir_sourceside =
