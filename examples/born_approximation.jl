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

# compute the sparse matrix correspinding to the FD stencil
ofds = ObsorbFDStencil(params);
rfds = RigidFDStencil(params);

# source function
src = Source(30, 1500, params; amp=100000, location_flag="distance");

# compute source side wavefield
dpdt= get_sourceside_wavefield(src, ofds, params);
imshow(dpdt[:,:,500]);
dpdt= reshape(dpdt, params.nz*params.nx, params.nt);

# forward born approximation
delta_lambda = zeros(params.data_format, params.nz*params.nx);

# put a single scatter in the middle of the model
iz = 50; ix=150; idx = (ix-1)*params.nz + iz;
delta_lambda[idx] = 1.0;

# pre-allocate recordings for scatter wavefield
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, params);

# forward modelling
born_approximation_forward!(rec, dpdt, delta_lambda, ofds, params);
gradient = born_approximation_adjoint(rec, dpdt, ofds, params);
gradient = reshape(gradient, params.nz, params.nx);

# dot_product test for born approximation
delta_lambda = randn(params.data_format, params.nz*params.nx);
born_approximation_forward!(rec, dpdt, delta_lambda, ofds, params);

# band-limited recordings
hl   = floor(Int64, length(src.p)/2)
rec1 = Recordings(irz, irx, params);
idx_o= [1]
for i = 1 : rec1.nr
    tmp = conv(randn(params.data_format, params.nt)*1000, src.p)
    copyto!(rec1.p, idx_o[1], tmp, hl+1, params.nt)
    idx_o[1] = idx_o[1] + params.nt
end
delta_lambda1 = born_approximation_adjoint(rec1, dpdt, ofds, params);

tmp1 = dot(delta_lambda, delta_lambda1)
tmp2 = dot(vec(rec.p), vec(rec1.p))
