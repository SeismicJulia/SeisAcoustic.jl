using SeisAcoustic

dir_work = joinpath(homedir(), "Desktop/Marmousi2");

path_vp   = joinpath(dir_work, "vp.rsf");
path_rho  = joinpath(dir_work, "rho.rsf");
path_wlet = joinpath(dir_work, "wavelet_08ms.rsf");

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
for i = 1 : ns
    l1 = (i-1)*delta_s
    l2 = l1 + delta_r*(nr-1)
    rec_loc = l1 : delta_r : l2
    irx[i] = floor.(Int64, rec_loc/dz) .+ 1
    irz[i] = 4 * ones(Int64, nr)
    isx_near[i] = floor(Int64, (l2+150     )/dz) + 1
    isz_near[i] = 3
    isx_far[i]  = floor(Int64, (l2+150+6000)/dz) + 1
    isz_far[i]  = 3
end

# near offset source
srcs    = get_multi_sources(isz_near, isx_near, fidiff; p=w);
dir_obs = joinpath(dir_work, "near_boat")
get_observations(dir_obs, irz, irx, srcs, fidiff)

# near offset source
srcs = get_multi_sources(isz_far, isx_far, fidiff; p=w);
dir_obs = joinpath(dir_work, "far_boat")
get_observations(dir_obs, irz, irx, srcs, fidiff)
