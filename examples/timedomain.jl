using SeisAcoustic

# velocity and density model
vel = 3000 * ones(100, 300);  # m/s
rho = 2000 * ones(100, 300);  # kg/m^3

# number of PML layers
npml = 20

# top boundary condition
free_surface = false    #(pml or free_surface)

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;  # use second as unit

# organize these parameters into a structure
params = ModelParams(rho, vel, npml, free_surface, dz, dx, dt, tmax;
         data_format=Float32, fd_flag="taylor", order=3)
# data_format = (Float32 or Float64)
# fd_flag     = ("taylor" or "ls")
# order       = 2 - 10 if we use "ls" to compute the FD coefficients
#             = 2 - n  if we use "taylor" expansion to compute FD coefficients

# compute the sparse matrix correspinding to the FD stencil
ofds = ObsorbFDStencil(params);

# check source related function
isz = 2; isx = 10;
src = Source(isz, isx, params)
src = Source(isz, isx, params;
      location_flag="distance", amp=4, type_flag="miniphase")
wavelet = rand(100)
src = Source(isz, isx, params; p=wavelet)

# multi-sources
isz = collect(1:3);
isx = collect(1:3);
srcs = get_multi_sources(isz, isx, params; amp=[1,2,3], type_flag=["ricker","miniphase", "sinc"])

# add source
spt = Snapshot(params)
add_source!(spt, src, 50)

wfd = Wavefield(params)
add_source!(wfd, src, 50)

spt = Snapshot(params)
add_multi_sources!(spt, srcs, 50)

(tmin, tmax) = time_range_multisources(srcs)

spt1 = Snapshot(params)
spt2 = Snapshot(params)
tmp1 = zeros(params.data_format, params.Nz * params.Nx)
tmp2 = zeros(params.data_format, params.Nz * params.Nx)

for it = 2 : 300

    one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
    add_source!(spt2, src, it)

    # prepare for next iteration
    copy_snapshot!(spt1, spt2)
end











# dominant frequency of source wavelet
fdom = 10;       #Hz
