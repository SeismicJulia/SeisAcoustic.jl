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
vel = vel[1:2:300,4000:2:4800];
vel = model_smooth(vel, 10);
rho = rho[1:2:300,4000:2:4800];

# SeisPlotTX(vel, hbox=2.0*2, wbox=5*2, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
# SeisPlotTX(rho, hbox=2.0*2, wbox=5*2, cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0008; tmax = 2.0;

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

# initialize a source
# isz = 2; isx = 50;
# src = Source(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

isx = collect(50:80:params.nx-50); isz = 2*ones(length(isx));
src = get_multi_sources(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# receiver location
irx = collect(1:2:params.nx);
irz = 2 * ones(Int64, length(irx));

# the folder save reflections
path_rec = joinpath(work_dir, "recordings")
argument_collection = (receiver_z=irz, receiver_x=irx, source=src,
                       fidiff_hete=params, fidiff_homo=params_homo, path_out=path_rec);

# do simulation
get_reflections(path_rec, irz, irx, src, params, params_homo)

SeisPlotTX(refl.p, dy=dt, cmap="gray", wbox=10, hbox=10);

# plotting
SeisPlotTX(hcat(dobs.p, direct_arrival.p, reflection.p), dy=dt, cmap="gray", wbox=30, hbox=10);
