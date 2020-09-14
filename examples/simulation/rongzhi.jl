using SeisPlot, SeisAcoustic

# ==============================================================================
#             read physical model
dir_work = joinpath(homedir(), "Desktop");

path_rho = joinpath(dir_work, "physical_model/rho.rsf");
path_vel = joinpath(dir_work, "physical_model/vel.rsf");

# read the physical model
println("Loading velocity and density model");
(hdr_rho, rho) = read_RSdata(path_rho);
(hdr_vel, vel) = read_RSdata(path_vel);

# cropped model for imaging
z_range = 1:3:300;
x_range = 4000:2:4600;
vel = vel[z_range, x_range];
rho = rho[z_range, x_range];

# grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0008; tmax = 1.25;

# top absorbing boundary to aviod surface-related multiple
free_surface = true;
data_format  = Float32;
order        = 2;

# tdparams for generating observations
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# source location and firing time
isx = collect(5:40:params.nx); ns=length(isx); isz = 2*ones(ns);
ot  = 0.5*rand(ns);
src = get_multi_sources(isz, isx, params; amp=100000, ot=ot, fdom=15);

# receiver location
irx = collect(1: 2 : params.nx);
irz = 2 * ones(Int64, length(irx));

# forward modeling of simultaneous sources
rec = multi_step_forward(src, params; rz=irz, rx=irx);
SeisPlotTX(rec.p)
