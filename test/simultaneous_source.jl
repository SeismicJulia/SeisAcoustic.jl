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
    SeisPlotTX(c2[:,i,:], hbox=5.6, wbox=4.4, pclip=95, cmap="gray",
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
    SeisPlotTX(tmp, hbox=5.6, wbox=4.4, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=[]); tight_layout();
    path_fig = join([dir_work "/figure/boat2_crg_clean_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();

    noise = tmp - c2[:,i,:]
    SeisPlotTX(noise, hbox=5.6, wbox=4.4, pclip=95, cmap="gray",
               dx=1, dy=0.001, xticks=[], yticks=[]); tight_layout();
    path_fig = join([dir_work "/figure/boat2_interference_" "$i" ".pdf"]);
    savefig(path_fig, dpi=100); close();
end

# ==============================================================================
#                           PGS data
# ==============================================================================
using PyPlot, SeisPlot, SeisAcoustic
dir_work = joinpath(homedir(), "Desktop/PGS/PGS_data_binary");

path_shots   = joinpath(dir_work, "data_Float64.rsf");
(hdr, shots) = read_RSdata(path_shots);

path_time    = joinpath(dir_work, "time_shift_Float64.rsf");
(hdr, tau)   = read_RSdata(path_time);

# get the source wavelet
channel_idx = 41;
SeisPlotTX(shots[:,channel_idx,:], hbox=10, wbox=8, pclip=90, cmap="gray",
           dx=1, dy=0.016, xticks=1:50:256, yticks=3:3:18, ticksize=20); tight_layout();

d = shots[1:100,51,85:145];
SeisPlotTX(d, hbox=10, wbox=8, pclip=99, cmap="gray",
           dx=1, dy=0.016); tight_layout();
wlet = sum(d, dims=2);
wlet = wlet[33:end];
taper= hanning(35*2+1); taper=taper[36:end]; wlet[33:end] .= wlet[33:end] .* taper;
figure(); plot(0.0:0.016:(length(wlet)-1)*0.016, wlet);
xlabel("Time (S)"); ylabel("Amplitude");
(fw, aw) = amplitude_spectra(wlet, 0.016); aw .= aw / maximum(aw);
(fd, ad) = amplitude_spectra(shots[:,:,1], 0.016); ad .= ad / maximum(ad);
figure(); plot(fw, aw, label="wavelet"); plot(fd, ad, label="shot gather");
xlabel("Frequency (Hz)"); ylabel("Relative amplitude");
legend();


# plot common shot gather
shot_idx = 27;
SeisPlotTX(shots[:,:,shot_idx], hbox=10, wbox=8, pclip=90, cmap="gray",
           dx=1, dy=0.016, xticks=1:50:256, yticks=3:3:18, ticksize=20); tight_layout();

# first arrival;
shot_idx = 75;
SeisPlotTX(shots[1:200,:,shot_idx], hbox=10, wbox=8, pclip=60, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=1:50:256, yticks=0.0:0.5:3.0, ticksize=20); tight_layout();

dx = 12.5
x = 150 : dx : 150+dx*255;
y = x / 1500
x1=collect(1:256);
plot(x1, y, "r");



d = shots[:,:,1];
d1    = trace_interpolation(d, 4);
d2    = trace_interpolation(d, 16);
SeisPlotTX(d, hbox=10, wbox=8, pclip=90, cmap="gray",
           dx=1, dy=0.016, xticks=1:50:256, yticks=3:3:18, ticksize=20); tight_layout();

SeisPlotTX(d1, hbox=10, wbox=8, pclip=90, cmap="gray",
           dx=1, dy=0.004, xticks=1:50:256, yticks=3:3:18, ticksize=20); tight_layout();

SeisPlotTX(d2, hbox=10, wbox=8, pclip=90, cmap="gray",
           dx=1, dy=0.001, xticks=1:50:256, yticks=3:3:18, ticksize=20); tight_layout();

function apply_time_shift(d::Matrix{Tv}, dt, tau) where {Tv<:AbstractFloat}

    (n1, n2) = size(d)
    fd       = fft(d, 1)

    dw = 1/dt / n1
    nw = floor(Int64, n1/2)+1



end

n1    = 100;
n2    = n1 * order;
dt    = 0.016; dt1 = dt / order;
x     = 0.0 : dt  : (n1-1)*dt;
x1    = 0.0 : dt1 : (n2-1)*dt1;
plot(x, d[1:n1]); plot(x1, d1[1:n2]);






# random time dithering for the second source used when artificailly blending
ot = 1.75 * rand(ns);
for i = 1 : ns
    ot[i] = round(Int64, ot[i]/dt) * dt
end


# source term
isz = 3; isx=2;
src = Source(isz, isx, fidiff; p=w);

# receiver
irx = collect(2:4:nx); irz = 4*ones(Int64, length(irx));
rec = Recordings(irz, irx, fidiff);
@time multi_step_forward!(rec, src, fidiff; print_flag=true);
SeisPlotTX(rec.p, wbox=10, hbox=10, cmap="gray", pclip=90);

#
SeisPlotTX(vp, wbox=6.8, hbox=1.4, cmap="rainbow", vmin=minimum(vp), vmax=maximum(vp)); tight_layout();

SeisPlotTX(vp, wbox=13.6, hbox=2.8, cmap="rainbow", vmin=minimum(vp), vmax=maximum(vp),
           dx=0.00125, dy=0.00125, xticks=2:2:16, yticks=0.5:0.5:3.0, ticksize=20,
           ylabel="Z (km)", labelsize=20); tight_layout();








#
