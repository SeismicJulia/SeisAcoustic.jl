using SeisPlot, SeisAcoustic, PyPlot

# ==================================kevin model ================================
# working directory
dir_work = joinpath(homedir(), "Desktop/simultaneous_source_data");

# physical model
path_vel  = joinpath(dir_work, "velocity_model/kevin_vel.rsf");
(hdr, vel)= read_RSdata(path_vel);
vel       = vel[1:2:300,301:700];
(nz, nx)  = size(vel);

# constant density model
rho       = ones(nz, nx);

# spatial grid size
dz = 8; dx = 8;

# time step and dominant frequency
dt = 1e-3; fdom = 12.;

# free surface top boundary condition
free_surface = true;

# simulate OBN acquisition, receiver location
irx = collect(1:2:nx); nr = length(irx);
irz = 5 * ones(Int64, length(irx));

# source location
isx = collect(2:2:20); ns = length(isx);
isz = 2 * ones(ns);

# ==============================================================================
#                    conventional survey, recording length
# ==============================================================================
# specify directory to save shot gathers
dir_con = joinpath(dir_work, "conventional_survey");
mkdir(dir_con);

# maximum simulation time for a single shot
tmax = 3.5;

# discretize the model using finite difference method
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax; data_format=Float32);

# generate a vector of sources
srcs= get_multi_sources(isz, isx, fidiff; fdom=fdom);

# parallel simulation over shots
get_shotgather(dir_con, irz, irx, srcs, fidiff);

# formulate a cube from separated shot gathers
cube  = zeros(fidiff.data_format, fidiff.nt, nr, ns);
for i = 1 : ns
    path_obs     = join([dir_con "/recordings_" "$i" ".bin"])
    rec          = read_recordings(path_obs)
    cube[:,:,i] .= rec.p
end

# write the cube to hard drive
hdr = RegularSampleHeader(cube);
path_cube = joinpath(dir_con, "cube.rsf");
write_RSdata(path_cube, hdr, cube);

# ==============================================================================
#                simultaneouse source simulation (two boat case)
# ==============================================================================
# specify directory to save simultaneouse survey
dir_sim = joinpath(dir_work, "simultaneouse_survey");
mkdir(dir_sim);

# recording length for a single shot
rec_length = 3.5;
nsample    = floor(Int64,rec_length/dt)+1

# simulate two boat firing simultaneouely
ns_half = floor(Int64, ns/2);

# the first boat is shooting regularly with time interval 3.5s
ot1     = collect(0.0 : rec_length : rec_length*(ns_half-1));

# the second boad is shooting randomly with respect to the first vessel.
ot2     = ot1 + 1.5 * rand(ns_half);

# save the shooting time (will be used for deblending process)
path_ot1 = joinpath(dir_sim, "shooting_time1.rsf");
hdr1 = RegularSampleHeader(ot1); write_RSdata(path_ot1, hdr1, ot1);

path_ot2 = joinpath(dir_sim, "shooting_time2.rsf");
hdr2 = RegularSampleHeader(ot2); write_RSdata(path_ot2, hdr2, ot2);

# maximum recording time (long continous recording)
tmax   = rec_length * ns_half + 2.0;

# finite difference grid
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax; data_format=Float32);
srcs   = get_multi_sources(isz, isx, fidiff; ot=vcat(ot1, ot2), fdom=fdom);

# continous recordings
path_rec = joinpath(dir_sim, "blended_shot.bin")
multi_step_forward!(irz, irx, srcs, fidiff, path_shot=path_rec);

# pseduo-deblending for boat 1 (align shots gather from boat 1)
rec   = read_recordings(path_rec);
cube1 = zeros(fidiff.data_format, nsample, nr, ns_half);
for i = 1 : ns_half
    tl = ot1[i];
    tu = ot1[i]+rec_length;
    i1 = floor(Int64, tl/dt)+1
    i2 = floor(Int64, tu/dt)+1
    cube1[:,:,i] .= rec.p[i1:i2, :]
end
path_boat1 = joinpath(dir_sim, "boat1_comb.rsf")
hdr        = RegularSampleHeader(cube1);
write_RSdata(path_boat1, hdr, cube1);

# combing for boat 2
cube2 = zeros(fidiff.data_format, nsample, nr, ns_half);
for i = 1 : ns_half
    tl = ot2[i];
    tu = ot2[i]+rec_length;
    i1 = floor(Int64, tl/dt)+1
    i2 = floor(Int64, tu/dt)+1
    cube2[:,:,i] .= rec.p[i1:i2, :]
end
path_boat2 = joinpath(dir_sim, "boat2_comb.rsf")
hdr        = RegularSampleHeader(cube2);
write_RSdata(path_boat2, hdr, cube2);

# ==============================================================================
#                            plotting
# ==============================================================================
dir_fig = joinpath(dir_work, "figure"); mkdir(dir_fig);

# plot velocity model
SeisPlotTX(vel, wbox=8, hbox=5, cmap="rainbow", vmax=maximum(vel), vmin=minimum(vel),
           ox=0, dx=0.00625, oy=0, dy=0.00625, xlabel="X (km)", ylabel="Z (km)");
path_fig = join([dir_fig "/velocity.pdf"]);
savefig(path_fig, dpi=100); close();

# plot conventional shot gathers
index_set = [1,3,6,9]
for i in index_set
    path_obs = join([dir_con "/recordings_" "$i" ".bin"]);
    rec      = read_recordings(path_obs);

    SeisPlotTX(rec.p, hbox=8, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=0.0:0.5:3.5, ticksize=25,
               ylabel="Time (S)", labelsize=25); tight_layout();
    path_fig = join([dir_fig "/recordings_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# plot boat 1
path_comb1 = joinpath(dir_sim, "boat1_comb.rsf");
(hdr, c1)  = read_RSdata(path_comb1);

# pseduo-deblended common shot gather
index_set = [1, 3, 5];
for i in index_set
    SeisPlotTX(c1[:,:,i], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=0.0:0.5:3.5); tight_layout();
    path_fig = join([dir_fig "/boat1_csg_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# pseduo-deblended common receiver gather
index_set = [60, 120, 180];
for i in index_set
    SeisPlotTX(c1[:,i,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=0.0:0.5:3.5); tight_layout();
    path_fig = join([dir_fig "/boat1_crg_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# plot boat 2
path_comb1 = joinpath(dir_sim, "boat2_comb.rsf");
(hdr, c2)  = read_RSdata(path_comb2);

# pseduo-deblended common shot gather
index_set = [1, 3, 5];
for i in index_set
    SeisPlotTX(c2[:,:,i], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=0.0:0.5:3.5); tight_layout();
    path_fig = join([dir_fig "/boat2_csg_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# pseduo-deblended common receiver gather
index_set = [60, 120, 180];
for i in index_set
    SeisPlotTX(c2[:,i,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=0.0:0.5:3.5); tight_layout();
    path_fig = join([dir_fig "/boat2_crg_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end


# ==============================================================================
#                          make animation
# ==============================================================================
# dir_work = joinpath(homedir(), "Desktop");
#
# path_rho = joinpath(dir_work, "physical_model/rho.rsf");
# path_vel = joinpath(dir_work, "physical_model/vel.rsf");
#
# # read the physical model
# (hdr_rho, rho0) = read_RSdata(path_rho);
# (hdr_vel, vel0) = read_RSdata(path_vel);
#
# # # cropped model for imaging
# vel = vel0[1:4:1350,150:4:2300];
# rho = rho0[1:4:1350,150:4:2300];
#
# # plotting velocity model
# SeisPlotTX(vel, wbox=8, hbox=5, cmap="rainbow", vmax=maximum(vel), vmin=minimum(vel),
#            ox=0, dx=0.00625, oy=0, dy=0.00625, xlabel="X (km)", ylabel="Z (km)");
# tight_layout(); savefig("/Users/wenlei/Desktop/vel.pdf"); close();
#
# # vertical and horizontal grid size
# dz = 6.25; dx = 6.25;
#
# # time step size and maximum modelling length
# dt = 0.0007; tmax = 6.0;
#
# # top boundary condition
# free_surface = true;
# data_format  = Float32;
# order        = 5;
#
# # tdparams for generating observations
# fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
#                   data_format=data_format, order=order);
#
# # receiver location
# irx = collect(1: 2: fidiff.nx-10);
# irz = 2 * ones(length(irx));
#
# # source location at isx=100
# isz = 2; isx = 90; ot=0.0;
# src = Source(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");
#
# # single source simulation
# interval = 10;  # save one every 10 snapshots
# path_pre = joinpath(dir_work, "single_pre.rsf")
# multi_step_forward!(src, fidiff, path_pre=path_pre, interval=interval);
#
# # single source recording
# rec = multi_step_forward!(irz, irx, src, fidiff);
# a = quantile(abs.(rec.p[:]), (98/100));
# SeisPlotTX(rec.p; wbox=4, hbox=5, cmap="seismic", vmax=a, vmin=-a,
#            ox=1, dx=1, oy=0.0, dy=0.0007, xlabel="Traces", ylabel="Time (s)")
# tight_layout(); savefig(joinpath(dir_work, "rec.pdf")); close();
#
#
# # wavefield animation
# (hdr, p) = read_RSdata(path_pre);
# a = quantile(abs.(p[:]), (98/100));
# path_ani = joinpath(dir_work, "wavefield_moive")
# SeisAnimation(path_ani, p; wbox=8, hbox=5, cmap="gray", vmax=a, vmin=-a,
#               ox=0, dx=0.00625, oy=0, dy=0.00625, title="wavefield", xlabel="X (km)", ylabel="Z (km)", interval=70)
#
# # make the record animation
# p = rec.p[1:interval:end,:];
# a = quantile(abs.(p[:]), (98/100));
# (nt, nr) = size(p);
# d = zeros(Float32, nt, nr, nt);
# for i = 1 : nt
#     d[1:i,:,i] .= p[1:i,:]
# end
# path_ani = joinpath(dir_work, "record_moive")
# SeisAnimation(path_ani, d; wbox=4, hbox=5, cmap="gray", vmax=a, vmin=-a,
#               ox=1, dx=1, oy=0.0, dy=0.007, xlabel="Traces", ylabel="Time (s)", interval=35)
#
#
# # simultaneouse source survey
# tmax = 12.0; # total simulation length
# fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
#                   data_format=data_format, order=order);
#
# # location of source
# isz = [2,2,2,2,2]; isx = [90, 180, 270, 360, 450];
#
# # start time
# ot=[0.0, 2.5, 3.1, 4.7, 5.3];
# srcs= get_multi_sources(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");
#
# interval = 20;
# path_pre = joinpath(dir_work, "multiple_pre.rsf")
# multi_step_forward!(srcs, fidiff, path_pre=path_pre, interval=interval);
#
# # single source recording
# rec = multi_step_forward!(irz, irx, srcs, fidiff);
#
#
# (hdr, p) = read_RSdata(path_pre);
# a = quantile(abs.(p[:]), (98/100));
# path_ani = joinpath(dir_work, "multiple_moive");
# SeisAnimation(path_ani, p; wbox=8, hbox=5, cmap="gray", vmax=a, vmin=-a,
#               ox=0, dx=0.00625, oy=0, dy=0.00625, title="wavefield", xlabel="X (km)", ylabel="Z (km)", interval=70)
#
# # make the record
# p = rec.p[1:interval:end,:];
# a = quantile(abs.(p[:]), (98/100));
# (nt, nr) = size(p);
# d = zeros(Float32, nt, nr, nt);
# for i = 1 : nt
#     d[1:i,:,i] .= p[1:i,:]
# end
# path_ani = joinpath(dir_work, "multiple_rec");
# SeisAnimation(path_rec, d; wbox=4, hbox=10, cmap="gray", vmax=a, vmin=-a,
#               ox=1, dx=1, oy=0.0, dy=0.014, xlabel="Traces", ylabel="Time (s)", interval=70)
