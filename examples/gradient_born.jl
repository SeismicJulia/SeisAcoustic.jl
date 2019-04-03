using SeisAcoustic, DSP, LinearAlgebra

# velocity and density model
vel = 3000 * ones(100, 300);  # m/s
rho = 2000 * ones(100, 300);  # kg/m^3
vel[51:100,:] .= 3100.;

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

# source function
src = Source(30, 1500, params; amp=100000, location_flag="distance");

# pre-allocate recordings for scatter wavefield
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, params);

# generate observations
multi_step_forward!(rec, src, ofds, params)

# initial velocity model
vel0  = 3000 * ones(100, 300);  # m/s
params0 = ModelParams(rho, vel0, npml, free_surface, dz, dx, dt, tmax;
         data_format=Float32, fd_flag="taylor", order=3);

# compute the sparse matrix correspinding to the FD stencil
ofds0 = ObsorbFDStencil(params0);

# generate synthetic
rec0 = Recordings(irz, irx, params);
multi_step_forward!(rec0, src, ofds0, params0)

# residue
res = Recordings(irz, irx, params);
for i = 1 : rec.nt * rec.nr
    res.p[i] = rec0.p[i] - rec.p[i]
end

# source-side wavefield
dpdt= get_sourceside_wavefield(src, ofds0, params0);
dpdt= reshape(dpdt, params.nz*params.nx, params.nt);

# the gradient
delta_lambda = born_approximation_adjoint(res, dpdt, ofds0, params0);
delta_lambda = reshape(delta_lambda, params.nz, params.nx)
