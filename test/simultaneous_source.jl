using SeisPlot, SeisAcoustic

dir_work = joinpath(homedir(), "Desktop");

path_rho = joinpath(dir_work, "physical_model/rho.rsf");
path_vel = joinpath(dir_work, "physical_model/vel.rsf");

# read the physical model
(hdr_rho, rho0) = read_RSdata(path_rho);
(hdr_vel, vel0) = read_RSdata(path_vel);

# # cropped model for imaging
vel = vel0[1:4:1350,150:4:2300];
rho = rho0[1:4:1350,150:4:2300];

SeisPlotTX(vel, wbox=8, hbox=5, cmap="rainbow", vmax=maximum(vel), vmin=minimum(vel)); tight_layout();
savefig("/Users/wenlei/Desktop/vel.pdf"); close();

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

# receiver location
irx = collect(1: 2: fidiff.nx-10);
irz = 2 * ones(length(irx));

# source location at isx=100
isz = 2; isx = 90; ot=0.0;
src = Source(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");

# single source simulation
path_pre = joinpath(dir_work, "single_pre.rsf")
multi_step_forward!(src, fidiff, path_pre=path_pre, interval=2);
rec = multi_step_forward!(irz, irx, src, fidiff);


function wave_animation(path_in::String, path_out="NULL";
         wbox=6, hbox=6, ox=0.0, dx=5.0, oz=0.0, dz=5.0,
         vmax="NULL", vmin="NULL", pclip=98,
         inter=5, cmap="PuOr", aspect="auto", interval=1)

    # get the header information
    hdr = read_RSheader(path)
    fid = open(path_pre, "r")

    # byte length of one pressure field
    pre_length = sizeof(hdr.data_format) * hdr.n1 * hdr.n2
    d          = zeros(hdr.data_format, hdr.n1, hdr.n2)

    fig = plt.figure(figsize=(wbox, hbox))
    ims = PyCall.PyObject[]

    for it = 1 : inter : nt

        # the location of pointer
        location = field_location[:data] + (it-1)*pre_length
        seek(fid, location)

        # read one pressure field
        read!(d, fid)

        # set range for plotting
        if vmax=="NULL" && vmin=="NULL"
           vmax = quantile(abs.(vec(d)), pclip/100.)
           vmin = -vmax
        end

        im = plt.imshow(d, cmap=cmap, vmin=vmin,vmax=vmax,extent=[ox, ox+(size(d,2)-1)*dx, oz+(size(d,1)-1)*dz,oz], aspect=aspect)
        push!(ims, PyCall.PyObject[im])
    end
    close(fid)
    ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
    pathout = join([pathout ".mp4"])
    ani[:save](pathout)
    return nothing
end


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


# ==============================================================================
#                       Gulf Mexico data
# ==============================================================================
scp wgao1@saig-ml.physics.ualberta.ca:/home/Data/GOM/DNCNN/data/train/input.bin /Users/wenlei/Desktop
scp wgao1@saig-ml.physics.ualberta.ca:/home/Data/GOM/DNCNN/data/train/output.bin /Users/wenlei/Desktop
scp wgao1@saig-ml.physics.ualberta.ca:/home/Data/GOM/DNCNN/data/train/pred.bin /Users/wenlei/Desktop


using MAT, PyPlot, SeisPlot, Statistics, SeisAcoustic

dir_radon = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Gulf_Mexico/Radon");
dir_ml    = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Gulf_Mexico/ML");
n1 = 1751; n2 = 183; n3 = 808;

path_delay = joinpath(dir_radon, "delay_time.mat");
dither     = matread(path_delay)["delay"];

dir_CCG    = joinpath(dir_radon, "CCG");
dir_CSG    = joinpath(dir_radon, "CSG");

# common shot gather (1751 * 183 * 200) [201-400]
path_csg_clean   = joinpath(dir_CSG, "D0.mat");
path_csg_pseudo  = joinpath(dir_CSG, "D1.mat");
path_csg_deblend = joinpath(dir_CSG, "D2.mat");

csg0 = matread(path_csg_clean)["D0"];
csg1 = matread(path_csg_pseudo)["D1"];
csg2 = matread(path_csg_deblend)["D2"];

# common channel gather (1751 * 4 * 808) [1, 61, 121, 181]
path_ccg_clean   = joinpath(dir_CCG, "D0.mat");
path_ccg_pseudo  = joinpath(dir_CCG, "D1.mat");
path_ccg_deblend = joinpath(dir_CCG, "D2.mat");

ccg0 = matread(path_ccg_clean)["D0"];
ccg1 = matread(path_ccg_pseudo)["D1"];
ccg2 = matread(path_ccg_deblend)["D2"];

# Machine learning
path_input = joinpath(dir_ml, "input.rsf");
(hdr, input) = read_RSdata(path_input);

path_output = joinpath(dir_ml, "output.rsf");
(hdr, output) = read_RSdata(path_output);

path_pred = joinpath(dir_ml, "pred.rsf");
(hdr, pred) = read_RSdata(path_pred);


# plotting common channel gather [1, 61, 121, 181]
cidx = 61; channel_idx = 2; sidx=1:2:n3;  v = quantile(abs.(vec(input[:, cidx, sidx])), 0.98);
h=10; w = 6; l = 25; t=25; xrange=80:80:400; yrange=1.5:1.5:6.0; dt=0.004; xl="Shot index"; yl="";
cmap="gray";
SeisPlotTX(input[:, cidx, sidx],                           hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg1_in.pdf", dpi=300); close();

SeisPlotTX(output[:, cidx, sidx],                          hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg1_true.pdf", dpi=300); close();

SeisPlotTX(pred[:, cidx, sidx],                            hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg1_pred.pdf", dpi=300); close();

SeisPlotTX(ccg2[:,channel_idx, sidx],                      hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg1_radon.pdf", dpi=300); close();

SeisPlotTX(output[:, cidx,sidx]-pred[:,cidx, sidx],        hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg1_pred_residue.pdf", dpi=300); close();

SeisPlotTX(ccg2[:,channel_idx, sidx]-output[:,cidx, sidx], hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg1_radon_residue.pdf", dpi=300); close();

# plotting common channel gather [1, 61, 121, 181]
cidx = 61; channel_idx = 2; sidx=2:2:n3;  v = quantile(abs.(vec(input[:, cidx, sidx])), 0.98);
h=10; w = 6; l = 25; t=25; xrange=80:80:400; yrange=1.5:1.5:6.0; dt=0.004; xl="Shot index"; yl="Time (s)";
cmap="gray";
SeisPlotTX(input[:, cidx, sidx],                           hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg2_in.pdf", dpi=300); close();

SeisPlotTX(output[:, cidx, sidx],                          hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg2_true.pdf", dpi=300); close();

SeisPlotTX(pred[:, cidx, sidx],                            hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg2_pred.pdf", dpi=300); close();

SeisPlotTX(ccg2[:,channel_idx, sidx],                      hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg2_radon.pdf", dpi=300); close();

SeisPlotTX(output[:, cidx,sidx]-pred[:,cidx, sidx],        hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg2_pred_residue.pdf", dpi=300); close();

SeisPlotTX(ccg2[:,channel_idx, sidx]-output[:,cidx, sidx], hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/ccg2_radon_residue.pdf", dpi=300); close();

# plotting common shot gather [201-400]
sidx=201; shot_idx =sidx-200; v = quantile(abs.(vec(input[:, :, sidx])), 0.99);
h=10; w = 6; l = 25; t=25; xrange=30:30:180; yrange=1.5:1.5:6.0; dt=0.004; xl="Channel index"; yl="";
cmap="gray";
SeisPlotTX(input[:, :, sidx],                       hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/csg_in.pdf", dpi=300); close();

SeisPlotTX(output[:, :, sidx],                      hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/csg_true.pdf", dpi=300); close();

SeisPlotTX(pred[:, :, sidx],                        hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/csg_pred.pdf", dpi=300); close();

SeisPlotTX(csg2[:, :, shot_idx],                    hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/csg_radon.pdf", dpi=300); close();

SeisPlotTX(output[:, :, sidx]-pred[:, :, sidx],     hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/csg_pred_residue.pdf", dpi=300); close();

SeisPlotTX(output[:, :, sidx]-csg2[:, :, shot_idx], hbox=h, wbox=w, vmin=-v, vmax=v, cmap=cmap, ylabel=yl, xlabel=xl,
           ox=1, dx=1, oy=0.0, dy=dt, xticks=xrange, yticks=yrange, ticksize=t, labelsize=l); tight_layout();
tight_layout(); savefig("/Users/wenlei/Desktop/csg_radon_residue.pdf", dpi=300); close();


# plot the cube of vessel 1
path="/Users/wenlei/Desktop/vessel1.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(input[:,:,1:2:n3]))); close(fid);
pscube < /Users/wenlei/Desktop/vessel1.bin size1=2.6 size2=2 size3=1.2 n1=1751 n2=183 n3=404 perc=0.98 > /Users/wenlei/Desktop/vessel1.eps

# SNR of DnCNN
r1 = compute_snr(output[:,:,201:400], pred[:,:,201:400]);
r2 = compute_snr(output[:,:,201:400], csg2);
r3 = compute_snr(output[:,:,201:400], input[:,:,201:400]);

function compute_snr(signal, result)

    ns = size(signal, 3)
    r  = zeros(ns)

    for i = 1 : ns
        r[i] = 10 * log10( norm(result[:,:,i]) / norm(signal[:,:,i]-result[:,:,i]))
    end

    return r
end

# ==============================================================================
#                       PGS data
# ==============================================================================
using PyPlot, SeisPlot, SeisAcoustic, FFTW

dir_pgs  = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/PGS/PGS_rsf");

path_data= joinpath(dir_pgs, "data.rsf");
(hdr, d) = read_RSdata(path_data);
(n1, n2, n3) = size(d);
dt = 0.016;


shot_idx = 100;
channel_idx = 80;

SeisPlotTX(d[:,:,shot_idx], hbox=10, wbox=8, pclip=95, cmap="gray", ylabel="Time (s)", xlabel="Channel",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=50:50:250, yticks=3:3:18, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/shot.pdf", dpi=300); close();

SeisPlotTX(d[:,channel_idx,:], hbox=10, wbox=8, pclip=95, cmap="gray", ylabel="Time (s)", xlabel="Channel",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=50:50:250, yticks=3:3:18, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/channel.pdf", dpi=300); close();

SeisPlotTX(d[233:357,channel_idx,80:120], hbox=5, wbox=3, pclip=95, cmap="gray", xticks=[], yticks=[]);
tight_layout(); savefig("/Users/wenlei/Desktop/room_channel.pdf", dpi=300); close();

# apply time shift to PGS data
path_data = joinpath(dir_pgs, "data_shifted.rsf");
(hdr, d_shift) = read_RSdata(path_data);

SeisPlotTX(d_shift[:,:,shot_idx], hbox=10, wbox=8, pclip=95, cmap="gray", ylabel="Time (s)", xlabel="Channel",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=50:50:250, yticks=3:3:18, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/shot.pdf", dpi=300); close();

SeisPlotTX(d_shift[:,channel_idx,:], hbox=10, wbox=8, pclip=95, cmap="gray", ylabel="Time (s)", xlabel="Channel",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=50:50:250, yticks=3:3:18, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/channel.pdf", dpi=300); close();

SeisPlotTX(d_shift[319:415,channel_idx,76:108], hbox=5, wbox=3, pclip=95, cmap="gray", xticks=[], yticks=[]);
tight_layout(); savefig("/Users/wenlei/Desktop/room_channel.pdf", dpi=300); close();

# estimate wavelet
SeisPlotTX(d[1:500,1,:], hbox=10, wbox=8, pclip=95, cmap="gray", ylabel="Time (s)", xlabel="Channel",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=50:50:250, yticks=2:2:8, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/ccg.pdf", dpi=300); close();

SeisPlotTX(d[1:100,1,90:145], hbox=10, wbox=8, pclip=95, cmap="gray",  ylabel="Time (s)", xlabel="Channel",
           ox=90, dx=1, oy=0.0, dy=0.016, xticks=100:10:140, yticks=0.3:0.3:1.5, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/cropped_wavelet.pdf", dpi=300); close();

w = sum(d[1:100,1,90:145], dims=2);
w1 = trace_interpolation(w, 4);
w2 = trace_interpolation(w1, 4);
SeisPlotTX(w1[29:end,:], hbox=10, wbox=4, xcur=0.9, style="wiggles", xlabel="Channel",
           xticks=[1], oy=0.0, dy=0.004, ox=1, dx=1, yticks=0.3:0.3:1.5, ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/wavelet.pdf", dpi=300); close();


(hdr, vel) = read_RSdata("/Users/wenlei/Desktop/lsrtm/physical_model/vel.rsf");
SeisPlotTX(vel, wbox=5.395*2, hbox=1.911*2, pclip=95, cmap="rainbow", vmax=maximum(vel), vmin=minimum(vel), xticks=[], yticks=[]);
tight_layout(); savefig("/Users/wenlei/Desktop/velocity.pdf", dpi=300); close();

(hdr, d) = read_RSdata(joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2/blend_4ms.rsf"));

SeisPlotTX(d[:,:,20], hbox=8, wbox=5.5, pclip=95, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/shot_clean.pdf", dpi=300); close();

SeisPlotTX(d[:,100,:], hbox=8, wbox=5.5, pclip=95, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/cog_clean.pdf", dpi=300); close();


# read pgs data
path_data= joinpath(dir_pgs, "data.rsf");
(hdr, pd) = read_RSdata(path_data);

# grad background noise
noise = d[1:140,end,:];
noise = trace_interpolation(noise, 4);

# generate synthetic noise
n1    = 5001;
(nn, n2) = size(noise);
n_new = zeros(n1, n2);
for i = 1 : size(noise,2)
    n_new[:,i] .= conv(noise[:,i], rand(n1).-0.5)[281:end-279];
end

# synthetic signal
(hdr, sd) = read_RSdata(joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2/blend_4ms.rsf"));
signal = sd[:,128,1:end-1];

(hdr, s1) = read_RSdata(joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2/near_4ms.rsf"));
output = s1[:,128,1:end-1];

# plotting
SeisPlotTX(pd[1:576,end,:], hbox=5, wbox=3, pclip=95, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.016, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/pgs_noise.pdf", dpi=300); close();

SeisPlotTX(noise, hbox=4, wbox=6, pclip=100, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/pure_noise.pdf", dpi=300); close();

SeisPlotTX(rand(n1,n2).-0.5, hbox=10, wbox=6, pclip=95, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/random_noise.pdf", dpi=300); close();

SeisPlotTX(n_new, hbox=10, wbox=6, pclip=99, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/synthetic_noise.pdf", dpi=300); close();


input = signal+n_new/5000;
output= output+n_new/5000;

SeisPlotTX(input[1:3000,:], hbox=10, wbox=6, pclip=96, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/input.pdf", dpi=300); close();

# SeisPlotTX(d[:,128,:], hbox=10, wbox=6, pclip=96, cmap="gray",
#            ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
# tight_layout(); savefig("/Users/wenlei/Desktop/check.pdf", dpi=300); close();

SeisPlotTX(output[1:3000,:], hbox=10, wbox=6, pclip=96, cmap="gray",
           ox=1, dx=1, oy=0.0, dy=0.004, xticks=[], yticks=[], ticksize=25, labelsize=25);
tight_layout(); savefig("/Users/wenlei/Desktop/output.pdf", dpi=300); close();


iend = 140;
SeisPlotTX(d[1:140,end,:], cmap="gray", pclip=95);





(hdr, d) = read_RSdata(joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2/blend_4ms.rsf"));
n1 = size(d, 1)









#
