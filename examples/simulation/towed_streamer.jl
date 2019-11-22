# ==============================================================================
#                                  BP model
# ==============================================================================
using SeisAcoustic

dir_work = joinpath(homedir(), "Desktop/BP");
dir_model= joinpath(dir_work, "physical_model");

dir_near = joinpath(dir_work, "near_boat");
dir_far  = joinpath(dir_work, "far_boat" );

path_vel  = joinpath(dir_model, "vel.rsf");
path_rho  = joinpath(dir_model, "rho.rsf");

(hdr, vel)  = read_RSdata(path_vel);
(hdr, rho)  = read_RSdata(path_rho);

# reduced physical model
vel = vel[1:1500, 700:end-700];
rho = rho[1:1500, 700:end-700];
(nz, nx) = size(vel);
for i2 = 1 : nx
    for i1 = 1 : nz
        if vel[i1,i2] > 4510
           vel[i1,i2] = 4510
        end
    end
end

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0008; tmax = 14.0;

# towed streamer
nr = 481; delta_r = 12.5; near_offset = 150.0;
ns = 510; delta_s = 25.0;

# compute receiver and source index for SLO geometry
irz = Vector{Vector{Int64}}(undef, ns);
irx = Vector{Vector{Int64}}(undef, ns);
isz_near = Vector{Int64}(undef, ns);
isx_near = Vector{Int64}(undef, ns);
isz_far  = Vector{Int64}(undef, ns);
isx_far  = Vector{Int64}(undef, ns);
for i = 1 : ns
    l2 = (nx-1) * dx - (i-1)*delta_s
    l1 = l2 - (nr-1) * delta_r
    rec_loc = l1 : delta_r : l2
    irx[i] = floor.(Int64, rec_loc/dx) .+ 1
    irz[i] = 4 * ones(Int64, nr)
    isx_near[i] = floor(Int64, (l1-150     )/dx) + 1
    isz_near[i] = 3
    isx_far[i]  = floor(Int64, (l1-150-6000)/dx) + 1
    isz_far[i]  = 3
end

dir_obs = joinpath(dir_work, "near_boat");
get_shotgather_adaptive(dir_obs, isz_near, isx_near, 1e5*ricker(10.0,dt), irz, irx, vel, rho, dz, dx, dt, tmax,
                        location_flag="index", data_format=Float64, order=2, free_surface=true, npad=100)

dir_obs = joinpath(dir_work, "far_boat");
get_shotgather_adaptive(dir_obs, isz_far, isx_far, 1e5*ricker(10.0,dt), irz, irx, vel, rho, dz, dx, dt, tmax,
                        location_flag="index", data_format=Float64, order=2, free_surface=true, npad=10)

scp -r wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/Desktop/BP/near_boat /Users/wenlei/Desktop/BP
rec = read_recordings("/Users/wenlei/Desktop/BP/near_boat/recordings_2.bin");
SeisPlotTX(rec.p, cmap="gray", wbox=10, hboz=10);


# download data to local machine
scp -r wgao1@hijitus.physics.ualberta.ca:/home/wgao1/Desktop/BP/far_boat $HOME/Desktop/BP/
scp -r wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/Desktop/BP/near_boat $HOME/Desktop/BP/

# resample the recordings from 0.8 ms to 4ms
dir_work = joinpath(homedir(), "Desktop/BP");
dt_new = 0.004; dt = 0.0008;
interval = floor(Int64, dt_new/dt)-1;

ns = 40;
dir_near = joinpath(dir_work, "near_boat");
dir_near_4ms = joinpath(dir_work, "near_boat_4ms");
for i =1 : ns
    path_in  = join([dir_near "/recordings_" "$i" ".bin"]);
    path_out = join([dir_near_4ms "/pressure_" "$i" ".rsf"]);
    rec      = read_recordings(path_in)
    p        = rec.p[1:interval:end,:]

    hdr      = RegularSampleHeader(p, d1=dt_new, d2=12.5)
    write_RSdata(path_out, hdr, p)
end

dir_far = joinpath(dir_work, "far_boat");
dir_far_4ms = joinpath(dir_work, "far_boat_4ms");
for i = 1 : ns
    path_in  = join([dir_far "/recordings_" "$i" ".bin"]);
    path_out = join([dir_far_4ms "/pressure_" "$i" ".rsf"]);
    rec      = read_recordings(path_in)
    p        = rec.p[1:interval:end,:]

    hdr      = RegularSampleHeader(p, d1=dt_new, d2=12.5)
    write_RSdata(path_out, hdr, p)
end

# make a cube
dir_near_4ms = joinpath(dir_work, "near_boat_4ms");
path = joinpath(dir_near_4ms, "pressure_1.rsf");
hdr  = read_RSheader(path);
p    = zeros(hdr.data_format, hdr.n1, hdr.n2, ns);
for i = 1 : ns
    path = join([dir_near_4ms "/pressure_" "$i" ".rsf"])
    (hdr, d) = read_RSdata(path)
    p[:,:,i].= d
end
SeisPlotTX(p[:,:,1], cmap="gray", wbox=10, hbox=10);

path = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2/near_4ms.rsf");
hdr  = RegularSampleHeader(p, d1=0.004, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
write_RSdata(path, hdr, p);

dir_far_4ms = joinpath(dir_work, "far_boat_4ms");
path = joinpath(dir_far_4ms, "pressure_1.rsf");
hdr  = read_RSheader(path);
p    = zeros(hdr.data_format, hdr.n1, hdr.n2, ns);
for i = 1 : ns
    path = join([dir_far_4ms "/pressure_" "$i" ".rsf"])
    (hdr, d) = read_RSdata(path)
    p[:,:,i].= d
end
path = joinpath(homedir(), "Desktop/Marmousi2/far_4ms.rsf");
hdr  = RegularSampleHeader(p, d1=0.004, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
write_RSdata(path, hdr, p);


# ==============================================================================
#                                   Marmousi2
# ==============================================================================
# using SeisAcoustic
#
# dir_work = joinpath(homedir(), "Desktop/Marmousi2");
# dir_model= joinpath(dir_work, "physical_model_rsf");
# dir_near = joinpath(dir_work, "near_boat");
# dir_far  = joinpath(dir_work, "far_boat" );
#
# path_vp   = joinpath(dir_model, "vp.rsf");
# path_rho  = joinpath(dir_model, "rho.rsf");
# path_wlet = joinpath(dir_model, "wavelet_08ms.rsf");
#
# (hdr, vp)   = read_RSdata(path_vp);
# (hdr, rho)  = read_RSdata(path_rho);
# (hdr, wlet) = read_RSdata(path_wlet);
#
# # reduced physical model
# vp  = vp[1:2:2401, 1000:5:end];
# rho = rho[1:2:2401, 1000:5:end];
# (nz, nx) = size(vp);
#
# # top boundary condition
# free_surface = true;
#
# # vertical and horizontal grid size
# dz = 6.25; dx = 6.25;
#
# # time step size and maximum modelling length
# dt = hdr.d1; tmax = 16.0;
#
# # precision
# data_format = Float32;
# order       = 2;
#
# # organize these parameters into a structure
# fidiff = TdParams(rho, vp, free_surface, dz, dx, dt, tmax;
#                   data_format=data_format, order=order);
#
# # towed streamer
# nr = 257; delta_r = 12.5; near_offset = 150.0;
# ns = 257; delta_s = 25.0;
#
# # compute receiver and source index for SLO geometry
# irz = Vector{Vector{Int64}}(undef, ns);
# irx = Vector{Vector{Int64}}(undef, ns);
# isz_near = Vector{Int64}(undef, ns);
# isx_near = Vector{Int64}(undef, ns);
# isz_far  = Vector{Int64}(undef, ns);
# isx_far  = Vector{Int64}(undef, ns);
# for i = 1 : ns
#     l2 = (nx-1) * dx - (i-1)*delta_s
#     l1 = l2 - (nr-1) * delta_r
#     rec_loc = l1 : delta_r : l2
#     irx[i] = floor.(Int64, rec_loc/dx) .+ 1
#     irz[i] = 4 * ones(Int64, nr)
#     isx_near[i] = floor(Int64, (l1-150     )/dx) + 1
#     isz_near[i] = 3
#     isx_far[i]  = floor(Int64, (l1-150-6000)/dx) + 1
#     isz_far[i]  = 3
# end
#
# # # # near offset source
# srcs    = get_multi_sources(isz_near, isx_near, fidiff; p=wlet);
# dir_obs = joinpath(dir_work, "near_boat");
# get_observations(dir_obs, irz, irx, srcs, fidiff);
#
# # far offset source
# il = 206; iu=ns;
# srcs = get_multi_sources(isz_far[il:iu], isx_far[il:iu], fidiff; p=wlet);
# dir_obs = joinpath(dir_work, "far_boat");
# get_observations(dir_obs, irz[il:iu], irx[il:iu], srcs, fidiff);
#
# # download data to local machine
# scp -r wgao1@neurus.physics.ualberta.ca:/$HOME/Desktop/Marmousi2/far_boat $HOME/Desktop/Marmousi2/
# scp -r wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/Desktop/Marmousi2/near_boat $HOME/Desktop/Marmousi2/
#
# # combine the file
# file_name = readdir("/Users/wenlei/Desktop/far_boat");
# dir_name  = joinpath(homedir(), "Desktop/rename")
# for i = 1 : length(file_name)
#     s = split(file_name[i], "_")[2]
#     a = parse(Int64, split(s, ".")[1]) + 205
#     path1 = joinpath("/Users/wenlei/Desktop/far_boat", file_name[i])
#     path2 = join([dir_name "/recordings_" "$a" ".bin"])
#     cp(path1, path2)
# end
#
# # resample the recordings from 0.8 ms to 4ms
# dir_work = joinpath(homedir(), "Desktop/Marmousi2");
# dt_new = 0.004; interval = floor(Int64, dt_new/dt)-1;
#
# dir_near = joinpath(dir_work, "near_boat");
# dir_near_4ms = joinpath(dir_work, "near_boat_4ms");
# for i = 1 : ns
#     path_in  = join([dir_near "/recordings_" "$i" ".bin"]);
#     path_out = join([dir_near_4ms "/pressure_" "$i" ".rsf"]);
#     rec      = read_recordings(path_in)
#     p        = rec.p[1:interval:end,:]
#
#     hdr      = RegularSampleHeader(p, d1=dt_new, d2=12.5)
#     write_RSdata(path_out, hdr, p)
# end
#
# dir_far = joinpath(dir_work, "far_boat");
# dir_far_4ms = joinpath(dir_work, "far_boat_4ms");
# for i = 1 : ns
#     path_in  = join([dir_far "/recordings_" "$i" ".bin"]);
#     path_out = join([dir_far_4ms "/pressure_" "$i" ".rsf"]);
#     rec      = read_recordings(path_in)
#     p        = rec.p[1:interval:end,:]
#
#     hdr      = RegularSampleHeader(p, d1=dt_new, d2=12.5)
#     write_RSdata(path_out, hdr, p)
# end
#
# # make a cube
# dir_work = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2")
# dir_near_4ms = joinpath(dir_work, "near_boat_4ms");
# path = joinpath(dir_near_4ms, "pressure_1.rsf");
# hdr  = read_RSheader(path);
# p    = zeros(hdr.data_format, hdr.n1, hdr.n2, ns);
# for i = 1 : ns
#     path = join([dir_near_4ms "/pressure_" "$i" ".rsf"])
#     (hdr, d) = read_RSdata(path)
#     p[:,:,i].= d
# end
# path = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2/near_4ms.rsf");
# hdr  = RegularSampleHeader(p, d1=0.004, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
# write_RSdata(path, hdr, p);
#
# dir_far_4ms = joinpath(dir_work, "far_boat_4ms");
# path = joinpath(dir_far_4ms, "pressure_1.rsf");
# hdr  = read_RSheader(path);
# p    = zeros(hdr.data_format, hdr.n1, hdr.n2, ns);
# for i = 1 : ns
#     path = join([dir_far_4ms "/pressure_" "$i" ".rsf"])
#     (hdr, d) = read_RSdata(path)
#     p[:,:,i].= d
# end
# path = joinpath(homedir(), "Desktop/Marmousi2/far_4ms.rsf");
# hdr  = RegularSampleHeader(p, d1=0.004, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
# write_RSdata(path, hdr, p);
#
#
# path1 = joinpath(homedir(), "Desktop/Marmousi2/near_4ms.rsf");
# path2 = joinpath(homedir(), "Desktop/Marmousi2/far_4ms.rsf");
# (hdr, p1) = read_RSdata(path1);
# (hdr, p2) = read_RSdata(path2);
#
# # manual blending the two cube
# delay_max = 1.75; tau = 1.75 * rand(ns);
# delay = zeros(Float64, ns);
# for i = 1 : ns
#     delay[i] = (round(Int64, tau[i]/dt)) * dt
# end
# hdr = RegularSampleHeader(delay);
# path = joinpath(homedir(), "Desktop/Marmousi2/delay_4ms.rsf");
# write_RSdata(path, hdr, delay);
#
# path = joinpath(homedir(), "Desktop/Marmousi2/delay_4ms.rsf");
# (hdr, delay) = read_RSdata(path);
# for i = 1 : ns
#     delay_idx = round(Int64, delay[i]/dt)
#     tmp = vcat(zeros(Float32, delay_idx, nr), p2[1:end-delay_idx,:,i])
#     p1[:,:,i] .= p1[:,:,i] .+ tmp
# end
# hdr = RegularSampleHeader(p1, d1=0.004, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
# path = joinpath(homedir(), "Desktop/Marmousi2/blend_4ms.rsf");
# write_RSdata(path, hdr, p1);
#
#
# # save it as binary
# path_in = joinpath(homedir(), "Desktop/Marmousi2/near_4ms.rsf");
# (hdr, near) = read_RSdata(path_in);
# path_out = joinpath(homedir(), "Desktop/Marmousi2/near_4ms.bin");
# fid  = open(path_out, "w");
# write(fid, vec(near)); close(fid);
#
# path_in = joinpath(homedir(), "Desktop/Marmousi2/far_4ms.rsf");
# (hdr, far) = read_RSdata(path_in);
# path_out = joinpath(homedir(), "Desktop/Marmousi2/far_4ms.bin");
# fid  = open(path_out, "w");
# write(fid, vec(far)); close(fid);
#
# path_in = joinpath(homedir(), "Desktop/Marmousi2/blend_4ms.rsf");
# (hdr, blend) = read_RSdata(path_in);
# path_out = joinpath(homedir(), "Desktop/Marmousi2/blend_4ms.bin");
# fid  = open(path_out, "w");
# write(fid, vec(blend)); close(fid);
#
# path_in = joinpath(homedir(), "Desktop/Marmousi2/delay_4ms.rsf");
# (hdr, dither) = read_RSdata(path_in);
# path_out = joinpath(homedir(), "Desktop/Marmousi2/dither_4ms.bin");
# fid  = open(path_out, "w");
# write(fid, convert(Vector{Float32}, dither)); close(fid);
#
#
# # make 16ms sampling cube
# p1_down = p1[1:4:end,:,:];
# p2_down = p2[1:4:end,:,:];
# path1 = joinpath(homedir(), "Desktop/Marmousi2/near_16ms.rsf");
# path2 = joinpath(homedir(), "Desktop/Marmousi2/far_16ms.rsf");
# hdr1  = RegularSampleHeader(p1_down, d1=0.016, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
# hdr2  = RegularSampleHeader(p2_down, d1=0.016, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
# write_RSdata(path1, hdr1, p1_down);
# write_RSdata(path2, hdr2, p2_down);
#
# shot_idx = 1
# SeisPlotTX(p1[:,:,shot_idx], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
# SeisPlotTX(p2[:,:,shot_idx], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
#
# channel_idx = 20
# SeisPlotTX(p1[:,channel_idx,:], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
# SeisPlotTX(p2[:,channel_idx,:], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
#
# # manual blending the two cube
# delay_max = 1.75; tau = 1.75 * rand(ns);
# delay = zeros(Float64, ns);
# for i = 1 : ns
#     delay[i] = (round(Int64, tau[i]/dt)) * dt
# end
# hdr = RegularSampleHeader(delay);
# path = joinpath(homedir(), "Desktop/Marmousi2/delay_16ms.rsf");
# write_RSdata(path, hdr, delay);
#
# path = joinpath(homedir(), "Desktop/Marmousi2/delay_16ms.rsf");
# (hdr, delay) = read_RSdata(path);
# for i = 1 : ns
#     delay_idx = round(Int64, delay[i]/dt)
#     tmp = vcat(zeros(Float32, delay_idx, nr), p2[1:end-delay_idx,:,i])
#     p1[:,:,i] .= p1[:,:,i] .+ tmp
# end
# hdr = RegularSampleHeader(p1, d1=0.016, d2=12.5, d3=25.0, unit1="s", unit2="meter", unit3="meter");
# path = joinpath(homedir(), "Desktop/Marmousi2/blend_16ms.rsf");
# write_RSdata(path, hdr, p1);
#
# for i = 60 : 60 : ns
#     SeisPlotTX(p1[:,:,i], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
# end
#
# for i = 40 : 40 : nr
#     SeisPlotTX(p1[:,i,:], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
# end
#
#
# # examine the data by plotting
# idx = 10;
# path_rec = join([dir_near_4ms "/pressure_" "$idx" ".rsf"]);
# (hdr, p) = read_RSdata(path_rec);
# SeisPlotTX(p, hbox=8, wbox=5, pclip=90, cmap="gray",
#            dx=1, dy=0.004, ox=1, oy=0.0, xticks=50:50:nr, yticks=2:2:16, ticksize=15,
#            xlabel="Trace number", ylabel="Time (S)", labelsize=15); tight_layout();
#
# # amplitude spectra
# (f, a) = amplitude_spectra(p, 0.004);
# plot(f, a);
#
# iter = 200;
# a    = 0.95;
# x    = zeros(iter);
# for i = 1 : iter
#     x[i] = a^i
# end
# plot(x)
#
# rec = read_recordings(joinpath(dir_obs, "recordings_8.bin"))
# SeisPlotTX(rec.p)
#
# # # # near offset source
# srcs    = get_multi_sources(isz_near[1:2], isx_near[1:2], fidiff; p=wlet);
# dir_obs = joinpath(dir_work, "near_boat");
# get_shotgather(dir_obs, irz[1:2], irx[1:2], srcs, fidiff);
#
# # near offset source
# srcs = get_multi_sources(isz_far, isx_far, fidiff; p=wlet);
# dir_obs = joinpath(dir_work, "far_boat");
# get_observations(dir_obs, irz, irx, srcs, fidiff);
