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

# # cropped model for imaging
z_range = 1:3:300;
x_range = 4000:2:4600;
vel = vel[z_range, x_range];
rho = rho[z_range, x_range];
vel = model_smooth(vel, 10);

# used for removing direct arrival
rho1 = copy(rho); rho1 .= minimum(rho);
vel1 = copy(vel); vel1 .= vel[1];

# SeisPlotTX(vel, hbox=1.5*2, wbox=3*2, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
# SeisPlotTX(rho, hbox=1.5*2, wbox=3*2, cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0008; tmax = 1.25;

# top absorbing boundary to aviod surface-related multiple
free_surface = false;
data_format  = Float32;
order        = 2;

# tdparams for generating observations
fidiff_hete = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                       data_format=data_format, order=order);

# tdparams for removing direct arrival
fidiff_homo = TdParams(rho1, vel1, free_surface, dz, dx, dt, tmax;
                       data_format=data_format, order=order);

# tdparams for imaging
fidiff = TdParams(rho1, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# vector of source
println("Set acquisition geometry");
isx = collect(30 : 70 : fidiff.nx-30); isz = 2*ones(length(isx));
src = get_multi_sources(isz, isx, fidiff; amp=100000, fdom=20, type_flag="miniphase");

# receiver location
irx = collect(1: 2 : fidiff.nx);
irz = 2 * ones(Int64, length(irx));

dir_obs        = joinpath(dir_work, "observations");
dir_sourceside = joinpath(dir_work, "sourceside");

# simulate observed data (direct arrival removed)
println("Simulate observations");
get_reflections(dir_obs, irz, irx, src, fidiff_hete, fidiff_homo);

println("Compute source-side wavefield and save its boundary");
get_wavefield_bound(dir_sourceside, src, fidiff);

# dobs = read_recordings(joinpath(dir_obs, "recordings_1.bin"));
# SeisPlotTX(dobs.p, dy=dt, cmap="gray", wbox=6, hbox=6);

# lsrtm
println("CGLS");
born_params = (irz=irz, irx=irx, dir_sourceside=dir_sourceside, fidiff=fidiff, normalization_flag=true, mute_index=10);
(x, his) = cgls(born_approximation, dir_obs; dir_work=dir_work, op_params=born_params, maxIter=25,
                d_axpby=recordings_axpby!, m_axpby=image_axpby!, d_norm=recordings_norm, m_norm=l2norm_rsf);

# check the result
println("plotting");
path = joinpath(dir_work, "sourceside/normalization.rsf");
(hdr, s) = read_RSdata(path);

path = joinpath(dir_work, "iterations/iteration_1.rsf");
(hdr, m) = read_RSdata(path);
p = m .* vec(s);
p = reshape(p, fidiff.nz, fidiff.nx);
p1= laplace_filter(p);
SeisPlotTX(p1, cmap="gray", wbox=10, hbox=5, title="RTM", name=joinpath(dir_work, "rtm.pdf"));

path = joinpath(dir_work, "iterations/iteration_10.rsf");
(hdr, m) = read_RSdata(path);
p = m .* vec(s);
p = reshape(p, fidiff.nz, fidiff.nx);
p1= laplace_filter(p);
SeisPlotTX(p1, cmap="gray", wbox=10, hbox=5, title="LSRTM_10", name=joinpath(dir_work, "lsrtm_10.pdf"));

path = joinpath(dir_work, "iterations/iteration_25.rsf");
(hdr, m) = read_RSdata(path);
p = m .* vec(s);
p = reshape(p, fidiff.nz, fidiff.nx);
p1= laplace_filter(p);
SeisPlotTX(p1, cmap="gray", wbox=10, hbox=5, title="LSRTM_25", name=joinpath(dir_work, "lsrtm_25.pdf"));

plot(his); savefig(joinpath(dir_work, "convergence.pdf")); close();
