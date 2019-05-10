using SeisPlot, SeisAcoustic

# velocity and density model
nz = 101; nx = 301
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3
rho[25:end,:] .= 2300;  # m/s

# top boundary condition
free_surface = true    #(pml or free_surface)

# grid size
h = 10

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;  # use second as unit

# organize these parameters into a structure
fparams = FdParams(rho, vel, h, free_surface; data_format=Float64);
tparams = ModelParams(rho, vel, free_surface, h, h, dt, tmax; data_format=Float64);

# single source
src = Source(2, 150, fparams; dt=dt, amp=100000);

# recordings
# initialize recordings
irx = collect(1 : 2 : fparams.nx);
irz = 2 * ones(length(irx));
rec_f = Recordings(irz, irx, dt, tmax, fparams);
rec_t = Recordings(irz, irx, tparams);

# generate recordings
get_recordings!(rec_f, src, fparams, print_flag=true);
multi_step_forward!(rec_t, src, tparams);

# forward modeling in frequency domain
pre_f = get_wavefield_FDFD(src, fparams; tmax=tmax, print_flag=true);

# compare against the time-domain modelling result
src = Source(2, 150, fparams; dt=0.001, amp=100000);
path  = joinpath(homedir(), "Desktop/boundary.rsb");
multi_step_forward!(path, src, tparams; save_flag="pressure");
(hdr, pre_t) = read_RSdata(path);

#plotting
it = 400
figure(figsize=(9,3)); imshow(pre_f[:,:,it*2], cmap="seismic");
figure(figsize=(9,3)); imshow(pre_t[:,:,it], cmap="seismic");

iz = 51; ix = 151;
figure(); plot(pre_f[iz,ix,:]); plot(pre_t[iz,ix,:])

# get recordigns in frequency domain
irx = collect(1:2:fparams.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, fparams);

# get Helmholtz operator
H = get_helmholtz_LU(fparams, 20.0);
b = zeros(Complex{Float64}, fparams.Nz * fparams.Nx);
idx = isx + fparams
