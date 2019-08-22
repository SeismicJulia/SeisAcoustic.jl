using SeisPlot, SeisAcoustic, LinearAlgebra

# homogeneous velocity and density model
nz = 101; nx = 301;
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# number of PML layers
npml = 20;

# top boundary condition
free_surface = false;   #(pml or free_surface)

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;  # use second as unit

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=Float32, fd_flag="taylor", order=10, npml=20, apml=900.);

# shows the default value for the keyword parameters
# data_format = (Float32 or Float64)
# fd_flag     = ("taylor" or "ls")
# order       = 2 - 10 if we use "ls" to compute the FD coefficients
#             = 2 - n  if we use "taylor" expansion to compute FD coefficients

# initialize a source
src = Source(35, 150, params; ot=0.0, fdom=20.0,
      type_flag="ricker", amp=100000, location_flag="index");

# initialize multi-sources
# isx = collect(5:60:295); ns=length(isx); isz = 2*ones(ns);
# ot  = 0.5*rand(ns);
# srcs = get_multi_sources(isz, isx, params; amp=100000, ot=ot, fdom=15);

# initialize recordings
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, params);

# save the boundary value and last wavefield
path_bnd = joinpath(homedir(), "Desktop/tmp_data/boundary_value.rsb");
path_wfd = joinpath(homedir(), "Desktop/tmp_data/last_wavefield.rsb");

# forward modeling of simultaneous sources
multi_step_forward!(rec, src, params; path_bnd=path_bnd, path_wfd=path_wfd);
SeisPlotTX(rec.p, pclip=98);

# save the pressure field
path_pre = joinpath(homedir(), "Desktop/tmp_data/pressure.rsb");

# multi_step_forward!(path_pre, src, params; save_flag="pressure");
multi_step_forward!(path_pre, src, params; save_flag="pressure");
(hdr, p0) = read_RSdata(path_pre);

# reconstruct the pressure field forward
p1 = pressure_reconstruct_forward(path_bnd, src, params);

# reconstruct the pressure field backward
p2 = pressure_reconstruct_backward(path_bnd, path_wfd, src, params);

# testing wavefield reconstruction forward or backward
it = 500;
SeisPlotTX(p0[:,:,it], pclip=98, cmap="seismic", wbox=6, hbox=2);
SeisPlotTX(p1[:,:,it], pclip=98, cmap="seismic", wbox=6, hbox=2);
SeisPlotTX(p2[:,:,it], pclip=98, cmap="seismic", wbox=6, hbox=2);
norm(p0-p1) / norm(p0)
norm(p0-p2) / norm(p0)
