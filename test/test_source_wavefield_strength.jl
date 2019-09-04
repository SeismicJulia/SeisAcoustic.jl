using SeisPlot, SeisAcoustic, LinearAlgebra

# working directory
work_dir = joinpath(homedir(), "Desktop/tmp_LSRTM");

# size of model
nz = 101; nx = 301;

# 2-layers velocity model
vel = 3000 * ones(nz, nx);  # m/s
# vel[51:end,:] .= 3500;

# constant density model
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = false;

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;

# precision
data_format=Float64;
order=3

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# initialize a source
# isz = 2; isx = 150;
# src = Source(isz, isx, params; amp=100000, fdom=15);

isx = collect(5:10:295); ns=length(isx);
isz = 2   * ones(ns);
ot  = 0.5 * rand(ns);
pol = rand(Bool, ns);
amp = 100000 * ones(ns);
for i = 1 : ns
    if !pol[i]
       amp[i] = - amp[i]
    end
end
srcs = get_multi_sources(isz, isx, params; amp=amp, ot=ot, fdom=5);

# generate observed data
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
dobs= Recordings(irz, irx, params);

# generate recordings and boundary value
path_bnd = joinpath(work_dir, "bnd.bin")
path_wfd = joinpath(work_dir, "wfd.bin")
path_sws = joinpath(work_dir, "sws.bin")

@time multi_step_forward!(dobs, srcs, params; path_bnd=path_bnd, path_wfd=path_wfd, path_sws=path_sws);
@time pre = get_sourceside_wavefield(srcs, params);

(hdr, s1) = read_RSdata(path_sws);
s2 = zeros(nz, nx)
for i2 = 1 : nx
    for i1 = 1 : nz
        s2[i1,i2] = dot(pre[i1,i2,:], pre[i1,i2,:])
    end
end
norm(s1-s2) / norm(s1)

SeisPlotTX(dobs.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)
