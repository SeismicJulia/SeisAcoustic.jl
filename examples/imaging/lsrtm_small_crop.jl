using SeisPlot, SeisAcoustic

# ==============================================================================
#             read physical model
dir_work = joinpath(homedir(), "Desktop/lsrtm");
initialize_lsrtm(dir_work);

path_rho = joinpath(dir_work, "physical_model/rho.rsf");
path_vel = joinpath(dir_work, "physical_model/vel.rsf");

# read the physical model
(hdr_rho, rho) = read_RSdata(path_rho);
(hdr_vel, vel) = read_RSdata(path_vel);

# # cropped model for imaging
z_range = 1:2:300;
x_range = 4000:2:4600;
vel = vel[z_range, x_range];
rho = rho[z_range, x_range];
vel = model_smooth(vel, 10);

rho1 = copy(rho); rho1 .= minimum(rho);
vel1 = copy(vel); vel1 .= vel[1];

SeisPlotTX(vel, hbox=1.5*2, wbox=3*2, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
SeisPlotTX(rho, hbox=1.5*2, wbox=3*2, cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0008; tmax = 1.5;

# top boundary condition
free_surface = false;
data_format  = Float32;
order        = 5;

# tdparams for generating observations
fidiff_hete = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# tdparams for removing direct arrival
fidiff_homo = TdParams(rho1, vel1, free_surface, dz, dx, dt, tmax;
                       data_format=data_format, order=order);

# tdparams for imaging
fidiff = TdParams(rho1, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# a single source
# isz = 2; isx = params.nx;
# src = Source(isz, isx, params; amp=100000, fdom=20, type_flag="miniphase");

# vector of source
# isx = collect(10 : 15 : fidiff.nx-10); isz = 2*ones(length(isx)); 
isx = collect(30 : 70 : fidiff.nx-30); isz = 2*ones(length(isx));
src = get_multi_sources(isz, isx, fidiff; amp=100000, fdom=20, type_flag="miniphase");

# receiver location
irx = collect(1: 2 : fidiff.nx);
irz = 2 * ones(Int64, length(irx));

dir_obs        = joinpath(dir_work, "observations");
dir_sourceside = joinpath(dir_work, "sourceside");

get_reflections(dir_obs, irz, irx, src, fidiff_hete, fidiff_homo);
get_wavefield_bound(dir_sourceside, src, fidiff);

# dobs = read_recordings(joinpath(dir_obs, "recordings_4.bin"));
# SeisPlotTX(dobs.p, dy=dt, cmap="gray", wbox=6, hbox=6);

# parameters for born approximation
born_params = (irz=irz, irx=irx, dir_sourceside=dir_sourceside, fidiff=fidiff);

# input and output for born approximation
path_m1 = joinpath(dir_work, "model_space/image1.rsf");
path_m2 = joinpath(dir_work, "model_space/image2.rsf");
dir_d1  = joinpath(dir_work, "data_space/born_forward");
dir_d2  = dir_obs;

image  = 1000.0 * randn(fidiff.data_format, fidiff.nz*fidiff.nx);
hdr    = RegularSampleHeader(image, title="random image");
write_RSdata(path_m1, hdr, image);

# apply operation
born_approximation(dir_d1, path_m1, 1; born_params...);
born_approximation(path_m2, dir_d2, 2; born_params...);

# test
(hdr, m1) = read_RSdata(path_m1);
(hdr, m2) = read_RSdata(path_m2);
tmp1 = dot(m1, m2);

tmp2 = [0.0]
for i = 1 : length(src)
    file_name = file_name = join(["recordings_" "$i" ".bin"])
    rec1 = read_recordings(joinpath(dir_d1, file_name))
    rec2 = read_recordings(joinpath(dir_d2, file_name))
    tmp2[1] += dot(rec1.p, rec2.p)
end
(tmp2[1]-tmp1) / tmp2[1]


# plot
# rec = read_recordings(joinpath(dir_d1, "recordings_1.bin"));
# SeisPlotTX(rec.p, dy=dt, cmap="gray", wbox=6, hbox=6);
# (hdr, m2) = read_RSdata(path_m2);
# SeisPlotTX(reshape(m2, fidiff.nz, fidiff.nx), cmap="gray", wbox=10, hbox=5);
# m3 = laplace_filter(reshape(m2, fidiff.nz, fidiff.nx));
# SeisPlotTX(m3, cmap="gray", wbox=10, hbox=5, pclip=96);
#
# tmp = zeros(fidiff.data_format, fidiff.nz, fidiff.nx)
# for i2 = 1 : fidiff.nx
#     for i1 = 2 : fidiff.nz
#         tmp[i1,i2] = (rho[i1,i2]-rho[i1-1,i2]) / (rho[i1,i2]+rho[i1-1,i2])
#     end
# end
# SeisPlotTX(tmp, cmap="gray", wbox=10, hbox=5);


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
