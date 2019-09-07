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
isz = 2; isx = 50;
src = Source(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# vector of source
# isx = collect(30:70:params.nx-30); isz = 2*ones(length(isx));
# src = get_multi_sources(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# receiver location
irx = collect(1:1:params.nx);
irz = 2 * ones(Int64, length(irx));

# do simulation
path_rec = joinpath(work_dir, "recordings")
get_reflections(path_rec, irz, irx, src, params, params_homo)

# plotting the data
refl = read_recordings(joinpath(path_rec, "reflections_4.bin"));
SeisPlotTX(refl.p, dy=dt, cmap="gray", wbox=6, hbox=6);


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

path_sourceside = joinpath(work_dir, "sourceside");
get_wavefield_bound(path_sourceside, src, params);

(hdr1, d1) = read_RSdata(joinpath(path_sourceside, "normalization.rsf"));
SeisPlotTX(d1, dy=dt, cmap="rainbow", wbox=6, hbox=3);

(h1, s1) = read_RSdata(joinpath(path_sourceside, "strength_1.rsf"));
(h2, s2) = read_RSdata(joinpath(path_sourceside, "strength_2.rsf"));
(h3, s3) = read_RSdata(joinpath(path_sourceside, "strength_3.rsf"));
(h4, s4) = read_RSdata(joinpath(path_sourceside, "strength_4.rsf"));
SeisPlotTX(s1, dy=dt, cmap="rainbow", wbox=6, hbox=3);
