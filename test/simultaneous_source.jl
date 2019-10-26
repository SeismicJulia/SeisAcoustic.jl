using SeisPlot, SeisAcoustic

dir_work = joinpath(homedir(), "Desktop/lsrtm");

path_rho = joinpath(dir_work, "physical_model/rho.rsf");
path_vel = joinpath(dir_work, "physical_model/vel.rsf");

# read the physical model
(hdr_rho, rho0) = read_RSdata(path_rho);
(hdr_vel, vel0) = read_RSdata(path_vel);

# # cropped model for imaging
vel = vel0[1:4:1350,650:4:2300];
rho = rho0[1:4:1350,650:4:2300];

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0007; tmax = 5.0;

# top boundary condition
free_surface = true;
data_format  = Float32;
order        = 5;

# tdparams for generating observations
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# initialize recordings
irx = collect(102:2:302);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, fidiff);

# source location at isx=100
isz = 2; isx = 90; ot=0.0;
src = Source(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");
path_pre = joinpath(dir_work, "wavefield/pressure_90.rsf");
multi_step_forward!(rec, src, fidiff, path_pre=path_pre, interval=20);
path_shot = joinpath(dir_work, "wavefield/shot_90.rec");
write_recordings(path_shot, rec);

# source location at isx=300
isx = 310;
src = Source(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");
path_pre = joinpath(dir_work, "wavefield/pressure_310.rsf");
multi_step_forward!(rec, src, fidiff, path_pre=path_pre, interval=20);
path_shot = joinpath(dir_work, "wavefield/shot_310.rec");
write_recordings(path_shot, rec);

# simultaneouse source survey
tmax = 7.0;
fidiff1 = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

isz = [2,2]; isx = [90, 310]; ot=[0.0, 2.0];
srcs= get_multi_sources(isz, isx, fidiff1; ot=ot, amp=100000, fdom=20, type_flag="miniphase");

# initialize recordings
irx = collect(102:2:302);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, fidiff1);

path_pre = joinpath(dir_work, "wavefield/pressure.rsf");
multi_step_forward!(rec, srcs, fidiff1, path_pre=path_pre, interval=20);
path_shot = joinpath(dir_work, "wavefield/shot.rec");
write_recordings(path_shot, rec);

# plotting
SeisPlotTX(vel, hbox=3.38*3, wbox=4.13*3, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel),
           dx=0.00625, dy=0.00625, xticks=0:0.5:2.5, yticks=0.0:0.5:2.0, ticksize=25,
           xlabel="X (km)", ylabel="Z (km)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/vel.pdf")
savefig(path_fig, dpi=100); close();

path_pre = joinpath(dir_work, "wavefield/pressure_90.rsf");
(hdr, d) = read_RSdata(path_pre);
SeisPlotTX(d[:,:,60], hbox=3.38*3, wbox=4.13*3, pclip=98, cmap="gray",
           dx=0.00625, dy=0.00625, xticks=0:0.5:2.5, yticks=0.0:0.5:2.0, ticksize=25,
           xlabel="X (km)", ylabel="Z (km)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/pressure_90.pdf")
savefig(path_fig, dpi=100); close();

path_rec = joinpath(dir_work, "wavefield/shot_90.rec");
rec      = read_recordings(path_rec);
SeisPlotTX(rec.p, hbox=3.38*2.14285, wbox=4.13*1.25, pclip=95, cmap="gray",
           dx=1, dy=0.0007, xticks=0:30:100, yticks=1:1:5, ticksize=25,
           xlabel="Trace number", ylabel="Time (S)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/shot_90.pdf")
savefig(path_fig, dpi=100); close();

path_pre = joinpath(dir_work, "wavefield/pressure_310.rsf");
(hdr, d) = read_RSdata(path_pre);
SeisPlotTX(d[:,:,60], hbox=3.38*3, wbox=4.13*3, pclip=98, cmap="gray",
           dx=0.00625, dy=0.00625, xticks=0:0.5:2.5, yticks=0.0:0.5:2.0, ticksize=25,
           xlabel="X (km)", ylabel="Z (km)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/pressure_310.pdf")
savefig(path_fig, dpi=100); close();

path_rec = joinpath(dir_work, "wavefield/shot_310.rec");
rec      = read_recordings(path_rec);
SeisPlotTX(rec.p, hbox=3.38*2.14285, wbox=4.13*1.25, pclip=95, cmap="gray",
           dx=1, dy=0.0007, xticks=0:30:100, yticks=1:1:5, ticksize=25,
           xlabel="Trace number", ylabel="Time (S)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/shot_310.pdf")
savefig(path_fig, dpi=100); close();

path_pre = joinpath(dir_work, "wavefield/pressure.rsf");
(hdr, d) = read_RSdata(path_pre);
SeisPlotTX(d[:,:,200], hbox=3.38*3, wbox=4.13*3, pclip=98, cmap="gray",
           dx=0.00625, dy=0.00625, xticks=0:0.5:2.5, yticks=0.0:0.5:2.0, ticksize=25,
           xlabel="X (km)", ylabel="Z (km)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/pressure.pdf")
savefig(path_fig, dpi=100); close();

path_rec = joinpath(dir_work, "wavefield/shot.rec");
rec      = read_recordings(path_rec);
SeisPlotTX(rec.p, hbox=3.38*3, wbox=4.13*1.25, pclip=95, cmap="gray",
           dx=1, dy=0.0007, xticks=0:30:100, yticks=1:1:7, ticksize=25,
           xlabel="Trace number", ylabel="Time (S)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "wavefield/shot.pdf")
savefig(path_fig, dpi=100); close();


# ==================================kevin model ================================
using PyPlot, SeisPlot, SeisAcoustic

# working directory
dir_work = joinpath(homedir(), "Desktop/simultaneous_source_data");

# physical model
path_vel = joinpath(dir_work, "velocity_model/kevin_vel.rsf");
(hdr,vel)= read_RSdata(path_vel); vel = vel[1:2:300,301:700]; (nz, nx) = size(vel);
rho      = ones(nz, nx);

# spatial and time grid size
dz = 8; dx = 8;
dt = 1e-3; fdom = 12.;

# top boundary condition
free_surface = true;

# simulate OBN acquisition, receiver location
irx = collect(1:2:nx); nr = length(irx);
irz = 2 * ones(length(irx));


# source location
isx = collect(2:2:nx); ns = length(isx);
isz = 2 * ones(ns);
srcs= get_multi_sources(isz, isx, fidiff; fdom=fdom);


# conventional survey, recording length
tmax = 3.5;
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax; data_format=Float32);
dir_obs = joinpath(dir_work, "conventional_survey");
get_observations(dir_obs, irz, irx, srcs, fidiff);

# make a cube from the conventional survey
cube  = zeros(fidiff.data_format, fidiff.nt, nr, ns);
for i = 1 : ns
    path_obs = join([dir_work "/conventional_survey/shot_" "$i" ".bin"])
    rec      = read_recordings(path_obs)
    cube[:,:,i] .= rec.p
end
hdr = RegularSampleHeader(cube);
path_cube = joinpath(dir_work, "conventional_survey/cube.rsf");
write_RSdata(path_cube, hdr, cube);


# shooting time for each shot
path_ot  = joinpath(dir_work, "simultaneous_source/fire_time.rsf");
(hdr, ot)= read_RSdata(path_ot);
tmax   = 352;
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax; data_format=Float32);
srcs   = get_multi_sources(isz, isx, fidiff; ot=ot, fdom=fdom);
rec    = Recordings(irz, irx, fidiff)

# continous recordings
path_rec = joinpath(dir_sim, "blended_shot.bin")
multi_step_forward!(rec, srcs, fidiff);
write_recordings(path_rec, rec);

# recording length for deblending
rec_length = 3.5;
ns_half    = floor(Int64, ns/2);

# combing for boat 1
c1 = zeros(fidiff.data_format, floor(Int64,rec_length/dt)+1, nr, ns_half);
for i = 1 : ns_half
    tl = ot[i]; tu = ot[i]+rec_length
    i1 = floor(Int64, tl/dt)+1
    i2 = floor(Int64, tu/dt)+1
    c1[:,:,i] .= rec.p[i1:i2, :]
end
path_boat1 = joinpath(dir_sim, "boat1_comb.rsf")
hdr        = RegularSampleHeader(c1);
write_RSdata(path_boat1, hdr, c1);

# combing for boat 2
c2 = zeros(fidiff.data_format, floor(Int64,rec_length/dt)+1, nr, ns_half);
for i = ns_half+1 : ns
    tl = ot[i]; tu = ot[i]+rec_length
    i1 = floor(Int64, tl/dt)+1
    i2 = floor(Int64, tu/dt)+1
    c2[:,:,i] .= rec.p[i1:i2, :]
end
path_boat2 = joinpath(dir_sim, "boat2_comb.rsf")
hdr        = RegularSampleHeader(c2);
write_RSdata(path_boat2, hdr, c2);


# ==============================================================================
#                            plotting
# ==============================================================================
# plot velocity model
# working directory
dir_work = joinpath(homedir(), "Desktop/simultaneous_source_data");

# physical model
path_vel = joinpath(dir_work, "velocity_model/kevin_vel.rsf");
(hdr,vel)= read_RSdata(path_vel); vel = vel[1:2:300,301:700]; (nz, nx) = size(vel);

SeisPlotTX(vel, hbox=2.5, wbox=21, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel),
           dx=0.008, dy=0.008, xticks=[], yticks=0.3:0.3:1.2, ticksize=25,
           ylabel="Z (km)", labelsize=25); tight_layout();
path_fig = joinpath(dir_work, "figure/vel.pdf")
savefig(path_fig, dpi=100); close();

# plot conventional shot gathers
index_set = [5, 10, 15, 20, 25, 105, 110, 115, 120, 125]
for i in index_set
    path_obs = join([dir_work "/conventional_survey/shot_" "$i" ".bin"]);
    rec      = read_recordings(path_obs);

    SeisPlotTX(rec.p, hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=1:1:3, ticksize=25,
               ylabel="Time (S)", labelsize=25); tight_layout();
    path_fig = join([dir_work "/figure/shot_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# plot boat 1 combed common receiver gather
path_comb1 = joinpath(dir_work, "simultaneous_source/boat1_comb.rsf");
(hdr, c1)  = read_RSdata(path_comb1);
index_set = [60, 120, 180];
for i in index_set
    SeisPlotTX(c1[:,i,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=[]); tight_layout();
    path_fig = join([dir_work "/figure/boat1_crg_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# plot boat 2 combed common receiver gather
path_comb2 = joinpath(dir_work, "simultaneous_source/boat2_comb.rsf");
(hdr, c2)  = read_RSdata(path_comb2);
index_set = [60, 120, 180];
for i in index_set
    SeisPlotTX(c2[:,i,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=[]); tight_layout();
    path_fig = join([dir_work "/figure/boat2_crg_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# plot boat 1 conventional common receiver gather
path_cube = joinpath(dir_work, "conventional_survey/cube.rsf");
(hdr, c) = read_RSdata(path_cube);
index_set = [60, 120, 180];
for i in index_set
    tmp = c[:,i,101:200]
    SeisPlotTX(tmp, hbox=5.5, wbox=4.4, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=[]); tight_layout();
    path_fig = join([dir_work "/figure/boat2_crg_clean_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();

    noise = tmp - c2[:,i,:]
    SeisPlotTX(noise, hbox=5.6, wbox=4.4, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=[]); tight_layout();
    path_fig = join([dir_work "/figure/boat2_interference_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end


# ns_half = floor(Int64, ns/2);
# d = zeros(Float32, fidiff.nt, ns_half);
# for i = 1 : 100
#     tmp = isx[i]
#     path_rec = join([dir_work "/shot_" "$tmp" ".rec"]);
#     rec = read_recordings(path_rec)
#     d[:,i] .= rec.p[:,60]
# end




SeisPlotTX(d, hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
           dx=1, dy=0.001, xticks=[], yticks=1:1:3, ticksize=25,
           ylabel="Time (S)", labelsize=25); tight_layout();
path_fig = join([dir_work "/clean_crg.pdf"]);
savefig(path_fig, dpi=100); close();








ns_half = floor(Int64, ns/2);
rec_length = 3.5
ot1 = collect(0.0 : rec_length : rec_length*(ns_half-1));

tau = 2.0 * rand(ns_half);
ot2 = ot1 + tau;
ot  = vcat(ot1, ot2);

# save firing time
hdr = RegularSampleHeader(ot);
path_ot = joinpath(dir_work, "fire_time.rsf");
write_RSdata(path_ot, hdr, ot);



irx = collect(1:2:nx); nr = length(irx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, fidiff);
multi_step_forward!(rec, srcs, fidiff; interval=10000);

path_rec = join([dir_work "/blended_shot" ".rec"]);
write_recordings(path_rec, rec);

for i = 1 : 4
    i1 = (i-1)*100*1000 + 1
    i2 = i1 + 10000
    ot = (i1-1)*0.001
    SeisPlotTX(rec.p[i1:i2,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
               dx=1, oy=(i1-1)*0.001, dy=0.001, xticks=[], yticks=ot:2:ot+10, ticksize=25,
               ylabel="Time (S)", labelsize=25); tight_layout();
    path_fig = join([dir_work "/blended_slice_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

i1 = 1
i2 = i1 + 15500
ot = (i1-1)*0.001
tu = (i2-1)*0.001
tmp = rec.p[i1:i2,:]; tmp = permutedims(tmp, [2,1]);
SeisPlotTX(tmp, wbox=3.38*3.5, hbox=3, pclip=95, cmap="gray",
           dy=1, ox=(i1-1)*0.001, dx=0.001, xticks=ot:3:tu, yticks=[], ticksize=25); tight_layout();
path_fig = join([dir_work "/blended_start.pdf"]);
savefig(path_fig, dpi=100); close();

i1 = 336501
i2 = i1 + 15500
ot = (i1-1)*0.001
tu = (i2-1)*0.001
tmp = rec.p[i1:i2,:]; tmp = permutedims(tmp, [2,1]);
SeisPlotTX(tmp, wbox=3.38*3.5, hbox=3, pclip=95, cmap="gray",
           dy=1, ox=(i1-1)*0.001, dx=0.001, xticks=ot:3:tu, yticks=[], ticksize=25); tight_layout();
path_fig = join([dir_work "/blended_end.pdf"]);
savefig(path_fig, dpi=100); close();



path="/Users/wenlei/Desktop/c1.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(c1))); close(fid);
pscube < /Users/wenlei/Desktop/c1.bin size1=6 size2=4 size3=2 labelsize=25 n1=3501 d1=0.001 d1num=1.0 f1num=1.0 label1="Time (s)" n2=200 d2num=50 f2num=50 label2="Receiver" n3=100 d3num=30 f3num=30 label3="shot" xbox=0.5 ybox=0.5 perc=0.98 > /Users/wenlei/Desktop/c1.eps

hdr = RegularSampleHeader(c1);
path_c1 = joinpath(dir_work, "boat1_comb.rsf");
write_RSdata(path_c1, hdr, c1);


SeisPlotTX(c1[:,j,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
           dx=1, dy=0.001, xticks=[], yticks=1:1:3, ticksize=25,
           ylabel="Time (S)", labelsize=25, title="Reciever 180", titlesize=25); tight_layout();
path_fig = join([dir_work "/crg_" "$j" ".pdf"]);
savefig(path_fig, dpi=100); close();

# comb for the boat2
c2 = zeros(fidiff.data_format, floor(Int64,rec_length/dt)+1, nr, ns_half);
for i = 1 : ns_half
    tl = ot[i+ns_half]; tu = ot[i+ns_half]+rec_length
    i1 = floor(Int64, tl/dt)+1
    i2 = floor(Int64, tu/dt)+1
    c2[:,:,i] .= rec.p[i1:i2, :]
end

SeisPlotTX(c2[:,j,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
           dx=1, dy=0.001, xticks=[], yticks=1:1:3, ticksize=25,
           ylabel="Time (S)", labelsize=25, title="Reciever 5", titlesize=25); tight_layout();
path_fig = join([dir_work "/crg2_" "$j" ".pdf"]);
savefig(path_fig, dpi=100); close();

path="/Users/wenlei/Desktop/c2.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(c2))); close(fid);
pscube < /Users/wenlei/Desktop/c2.bin size1=6 size2=4 size3=2 labelsize=25 n1=3501 d1=0.001 d1num=1.0 f1num=1.0 label1="Time (s)" n2=200 d2num=50 f2num=50 label2="Receiver" n3=100 d3num=30 f3num=30 label3="shot" xbox=0.5 ybox=0.5 perc=0.98 > /Users/wenlei/Desktop/c2.eps


hdr = RegularSampleHeader(c2);
path_c2 = joinpath(dir_work, "boat2_comb.rsf");
write_RSdata(path_c2, hdr, c2);

j = 100
SeisPlotTX(c2[:,j,:], hbox=3.38*2.5, wbox=5, pclip=95, cmap="gray",
           dx=1, dy=0.001, xticks=[], yticks=1:1:3, ticksize=25,
           ylabel="Time (S)", labelsize=25); tight_layout();


# make animation
# (hdr, d) = read_RSdata(path_wfd);
# path_ani = joinpath(dir_work, "wavefield/animation_100");
# SeisAnimation(path_ani, d; wbox=4.13, hbox=3.38, cmap="gray", vmin=-10, vmax=10, xticks=[], yticks=[], interval=20);
