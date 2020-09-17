using SeisPlot, PyPlot, SeisAcoustic

# create a simple two layers velocity model with constant density,
nz = 100; nx = 301;
vel = 2800 * ones(nz, nx);  # m/s
vel[31:end,:] .= 3500;      # m/s
rho = 2.0 * ones(nz, nx);   # kg/m^3
rho[61:end,:] .= 2.5;       # m/s

# plotting
SeisPlotTX(vel, wbox=9, hbox=5, title="velocity", cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
SeisPlotTX(rho, wbox=9, hbox=5, title="density", cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));

# number of PML layers
npml = 20;

# top boundary condition
free_surface = true;        #(pml or free_surface)

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 1.25;     # use second as unit

# organize these parameters into a structure
# shows the default value for the keyword parameters
# data_format = (Float32 or Float64)
# fd_flag     = ("taylor" or "ls")
# order       = 2 - 10 if we use "ls" to compute the FD coefficients
#             = 2 - n  if we use "taylor" expansion to compute FD coefficients
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=3, npml=20, apml=900.);
# type "?TdParams" to get more information about this structure

# ==============================================================================
#                          single source simulation
# ==============================================================================
isz = 2; isx = floor(Int64, params.nx/2);
# type_flag = "ricker" or "miniphase" or a vector
# location_flag = "index" or "distance"
src = Source(isz, isx, params; ot=0.0, fdom=15.0,
      type_flag="miniphase", amp=100000, location_flag="index");

# receiver location
irx = collect(1:2:params.nx);
irz = 2 * ones(Int64, length(irx));

# get one recordings
rec = multi_step_forward(src, params; rz=irz, rx=irx);
SeisPlotTX(rec.p, dy=0.001, title="one shot")
# type "multi_step_forward" to get more information about forward modelling

# save wavefield every 10 steps
# boundary of wavefield which are used for imaging
# snapshot(vz, vx, pz, px), strength of source side wavefield and so on
path_pre = joinpath(homedir(), "Desktop/pressure.rsf");
path_bnd = joinpath(homedir(), "Desktop/boundary.rsf");
multi_step_forward(src, params; path_pre=path_pre, path_bnd=path_bnd, interval=10);

# read the wavefield
(hdr, pre) = read_RSdata(path_pre);
path_anim  = joinpath(homedir(), "Desktop/wave_propagation");
vmin = minimum(pre[:,:,20]); vmax = maximum(pre[:,:,20]);
SeisAnimation(path_anim, pre; wbox=9, hbox=5, vmin=vmin, vmax=vmax, dx=10, dy=10)

# save every thing generated during forward modelling
# path_shot= joinpath(homedir(), "Desktop/recordings.bin" );
# path_spt = joinpath(homedir(), "Desktop/snapshot.rsf" );
# path_wfd = joinpath(homedir(), "Desktop/wavefield.rsf");
# path_pre = joinpath(homedir(), "Desktop/pressure.rsf" );
# path_bnd = joinpath(homedir(), "Desktop/boundary.rsf" );
# path_lwfd= joinpath(homedir(), "Desktop/last_wfd.rsf" );
# path_sws = joinpath(homedir(), "Desktop/strength.rsf" );
# multi_step_forward(src, params; rz=irz, rx=irx, path_shot=path_shot, path_spt=path_spt, path_wfd=path_wfd, path_pre=path_pre, interval=1,
#                    path_bnd=path_bnd, path_lwfd=path_lwfd, path_sws=path_sws);

# ==============================================================================
#                          simultaneouse source
# ==============================================================================
isx = collect(5:40:params.nx); ns=length(isx); isz = 2*ones(ns);
ot  = 0.3*rand(ns);
# ricker wavelet as the default one
src = get_multi_sources(isz, isx, params; amp=100000, ot=ot, fdom=15);

# forward modeling of simultaneous sources
path_pre = joinpath(homedir(), "Desktop/pressure_semi.rsf");
rec1 = multi_step_forward(src, params;rz=irz, rx=irx, path_pre=path_pre, interval=10);
SeisPlotTX(rec1.p, dy=0.001, title="semitaneous source")

(hdr, pre) = read_RSdata(path_pre);
path_anim  = joinpath(homedir(), "Desktop/wave_propagation_semi");
vmin = minimum(pre[:,:,50]); vmax = maximum(pre[:,:,50]);
SeisAnimation(path_anim, pre; wbox=9, hbox=5, vmin=vmin, vmax=vmax, dx=10, dy=10)

# ==============================================================================
#                multiple shot simulation
# ==============================================================================
isx = collect(120:30:params.nx); ns=length(isx); isz = 2*ones(ns);
src = get_multi_sources(isz, isx, params; amp=100000, fdom=15);

# simulate OBN acquisition with fixed receiver location
irx = collect(1:2:params.nx);
irz = 2 * ones(Int64, length(irx));

# do parallel simulation
dir_obs = joinpath(homedir(), "Desktop/shotgather")
get_shotgather(dir_obs, irz, irx, src, params);

path_shot = joinpath(dir_obs, "recordings_3.bin");
rec       = read_recordings(path_shot)
SeisPlotTX(rec.p, dy=0.001, title="OBN shot");

# simulate towed streamer acquisition with moving receiver
ns  = length(src)
irx = Vector{Vector{Int64}}(undef, ns)
irz = Vector{Vector{Int64}}(undef, ns)
for i = 1 : ns
    irx[i]= collect(isx[i]-102: isx[i]-2);
    irz[i]= 2 * ones(Int64, length(irx[i]));
end

# do parallel simulation
dir_obs = joinpath(homedir(), "Desktop/shotgather_tow")
get_shotgather(dir_obs, irz, irx, src, params);

path_shot = joinpath(dir_obs, "recordings_7.bin");
rec       = read_recordings(path_shot)
SeisPlotTX(rec.p, dy=0.001, title="towed streamer");

# end of documentation
