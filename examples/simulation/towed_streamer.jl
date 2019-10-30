using SeisAcoustic

dir_work = joinpath(homedir(), "Desktop/Marmousi2");
dir_model= joinpath(dir_work, "physical_model_rsf");
dir_near = joinpath(dir_work, "near_boat");
dir_far  = joinpath(dir_work, "far_boat" );

path_vp   = joinpath(dir_model, "vp.rsf");
path_rho  = joinpath(dir_model, "rho.rsf");
path_wlet = joinpath(dir_model, "wavelet_08ms.rsf");

(hdr, vp)   = read_RSdata(path_vp);
(hdr, rho)  = read_RSdata(path_rho);
(hdr, wlet) = read_RSdata(path_wlet);

# reduced physical model
vp  = vp[1:2:2401, 1000:5:end];
rho = rho[1:2:2401, 1000:5:end];
(nz, nx) = size(vp);

# top boundary condition
free_surface = true;

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = hdr.d1; tmax = 16.0;

# precision
data_format = Float32;
order       = 2;

# organize these parameters into a structure
fidiff = TdParams(rho, vp, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# towed streamer
nr = 257; delta_r = 12.5; near_offset = 150.0;
ns = 257; delta_s = 25.0;

# compute receiver and source index for SLO geometry
irz = Vector{Vector{Int64}}(undef, ns);
irx = Vector{Vector{Int64}}(undef, ns);
isz_near = Vector{Int64}(undef, ns);
isx_near = Vector{Int64}(undef, ns);
isz_far  = Vector{Int64}(undef, ns);
isx_far  = Vector{Int64}(undef, ns);
for i = 128 : ns
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

# # # near offset source
srcs    = get_multi_sources(isz_near, isx_near, fidiff; p=wlet);
dir_obs = joinpath(dir_work, "near_boat");
get_observations(dir_obs, irz, irx, srcs, fidiff);

# near offset source
srcs = get_multi_sources(isz_far, isx_far, fidiff; p=wlet);
dir_obs = joinpath(dir_work, "far_boat");
get_observations(dir_obs, irz, irx, srcs, fidiff);

# download data to local machine
scp -r wgao1@neurus.physics.ualberta.ca:/$HOME/Desktop/Marmousi2/far_boat /$HOME/Desktop/Marmousi2/
scp -r wgao1@saig-ml.physics.ualberta.ca:/$HOME/Desktop/Marmousi2/near_boat /$HOME/Desktop/Marmousi2/


# resample the recordings from 0.8 ms to 4ms
dt_new = 0.016; interval = floor(Int64, dt_new/dt)-1;
dir_near_4ms = joinpath(dir_work, "near_boat_4ms");
for i = 1 : 76
    path_in  = join([dir_near "/recordings_" "$i" ".bin"]);
    path_out = join([dir_near_4ms "/pressure_" "$i" ".rsf"]);
    rec      = read_recordings(path_in)
    p        = rec.p[1:interval:end,:]

    hdr      = RegularSampleHeader(p, d1=dt_new, d2=12.5)
    write_RSdata(path_out, hdr, p)
end


# get common channel gather
rec_idx= 125; src_max = 76;
path_shot = join([dir_near_4ms "/pressure_1.rsf"]);
hdr       = read_RSheader(path_shot);
ccp       = zeros(hdr.data_format, hdr.n1, src_max);
for i = 1 : src_max
    path_shot = join([dir_near_4ms "/pressure_" "$i" ".rsf"]);
    (hdr, p)  = read_RSdata(path_shot)
    ccp[:,i] .= p[:,rec_idx]
end
SeisPlotTX(ccp, hbox=8, wbox=5, pclip=90, cmap="gray",
           dx=1, dy=0.016, ox=1, oy=0.0, xticks=20:20:src_max, yticks=2:2:16, ticksize=15,
           xlabel="Trace number", ylabel="Time (S)", labelsize=15); tight_layout();

# pseduo-blending
delay_max = 1.75; tau = 1.75 * rand(ns);
delay_idx = zeros(Int64, ns);
for i = 1 : ns
    tau[i] = round(Int64, tau[i]/dt) + 1
end



# examine the data by plotting
idx = 10;
path_rec = join([dir_near_4ms "/pressure_" "$idx" ".rsf"]);
(hdr, p) = read_RSdata(path_rec);
SeisPlotTX(p, hbox=8, wbox=5, pclip=90, cmap="gray",
           dx=1, dy=0.004, ox=1, oy=0.0, xticks=50:50:nr, yticks=2:2:16, ticksize=15,
           xlabel="Trace number", ylabel="Time (S)", labelsize=15); tight_layout();

# amplitude spectra
(f, a) = amplitude_spectra(p, 0.004);
plot(f, a);
