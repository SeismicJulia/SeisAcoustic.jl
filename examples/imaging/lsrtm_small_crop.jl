using SeisPlot, SeisAcoustic

# ==============================================================================
#             read physical model
work_dir = joinpath(homedir(), "Desktop/lsrtm");

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

# tdparams for generating observations
params_hete = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# tdparams for removing direct arrival
rho1 = copy(rho); rho1 .= minimum(rho);
vel1 = copy(vel); vel1 .= vel[1];
params_homo = TdParams(rho1, vel1, free_surface, dz, dx, dt, tmax;
                       data_format=data_format, order=order);

# tdparams for imaging
params = TdParams(rho1, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# a single source
# isz = 2; isx = params.nx;
# src = Source(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# vector of source
isx = collect(30:70:params.nx-30); isz = 2*ones(length(isx));
src = get_multi_sources(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# receiver location
irx = collect(1:2:params.nx);
irz = 2 * ones(Int64, length(irx));

dir_obs        = joinpath(work_dir, "data_space/observation");
dir_sourceside = joinpath(work_dir, "sourceside");
dir_born       = joinpath(work_dir, "data_space/born_forward");
path_m         = joinpath(work_dir, "model_space/m.rsf");

# ==============================================================================
#                get reflections (direct arrival been removed)
# ==============================================================================
get_reflections(dir_obs, irz, irx, src, params_hete, params_homo);

# ==============================================================================
#               the model used for imaging
# ==============================================================================
get_wavefield_bound(dir_sourceside, src, params);

# ==============================================================================
#                  the adjoint of born approximation
# ==============================================================================
get_born_adjoint(path_m, dir_obs, dir_sourceside, params; remove_flag=false);

# ==============================================================================
#                  the forward of born approximation
# ==============================================================================
m = 1000.0 * randn(params.data_format, params.nz*params.nx);
get_born_forward(dir_born, irz, irx, m, dir_sourceside, params);

dobs = read_recordings(joinpath(dir_born, "recording_1.bin"));
SeisPlotTX(dobs.p, dy=dt, cmap="gray", wbox=6, hbox=6);

(hdr, m1) = read_RSdata(path_m);
tmp1 = dot(m, m1)

ns = length(src)
tmp2 = [0.0]
for i = 1 : ns
    file_name = file_name = join(["recording_" "$i" ".bin"])
    rec1 = read_recordings(joinpath(dir_born, file_name))
    rec2 = read_recordings(joinpath(dir_obs, file_name))
    tmp2[1] += dot(rec1.p, rec2.p)
end




(hdr, m) = read_RSdata(path_m);

SeisPlotTX(reshape(m, params.nz, params.nx), cmap="gray", wbox=6, hbox=3);

tmp = zeros(params.data_format, params.nz, params.nx)
for i2 = 1 : params.nx
    for i1 = 2 : params.nz
        tmp[i1,i2] = (rho[i1,i2]-rho[i1-1,i2]) / (rho[i1,i2]+rho[i1-1,i2])
    end
end
SeisPlotTX(tmp, cmap="gray", wbox=6, hbox=3);

p = laplace_filter(reshape(m, params.nz, params.nx));
SeisPlotTX(p , cmap="gray", wbox=6, hbox=3);
SeisPlotTX(p1, cmap="gray", wbox=6, hbox=3);

# plotting the data
# dobs = read_recordings(joinpath(dir_obs, "reflection_1.bin"));
# SeisPlotTX(dobs.p, dy=dt, cmap="gray", wbox=6, hbox=6);


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
