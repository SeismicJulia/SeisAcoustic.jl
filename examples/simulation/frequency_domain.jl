using PyPlot, SeisAcoustic

# velocity and density model
vel = 3000 * ones(101, 301);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(101, 301);  # kg/m^3

# number of PML layers
npml = 20

# top boundary condition
free_surface = true    #(pml or free_surface)

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;  # use second as unit

# organize these parameters into a structure
params = ModelParams(rho, vel, free_surface, dz, dx, dt, tmax; data_format=Float64)

# single source
src = Source(2, 150, params; amp=100000);

# forward modeling in time domain
path_pre = joinpath(homedir(), "Desktop/pressure.rsb");
multi_step_forward(path_pre, src, ofds, params);
(hdr, pre0) = read_RSdata(path_pre);

# forward modeling in frequency domain
pre1 = get_wavefield_FDFD(src, params; print_flag=true)

#plotting
it = 400
figure(); imshow(pre0[:,:,it]);
figure(); imshow(pre1[:,:,it]);
