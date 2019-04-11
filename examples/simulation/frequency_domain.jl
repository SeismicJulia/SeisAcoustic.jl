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
params = ModelParams(rho, vel, npml, free_surface, dz, dx, dt, tmax; data_format=Float32)

# frequency
f = 20.0;
omega = 2 * pi * f;

H = get_helmholtz_LU(params, omega);
