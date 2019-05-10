using SeisPlot, SeisAcoustic

# velocity and density model
nz = 101; nx = 301
vel = 3000 * ones(nz, nx);  # m/s
# vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = true    #(pml or free_surface)

# grid size
h = 10

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;  # use second as unit

# organize these parameters into a structure
params = FdParams(rho, vel, h, free_surface; data_format=Float64);

# single source
src = Source(2, 150, params; dt=dt, amp=100000);

# recordings
irx = collect(1 : 2 : fparams.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, dt, tmax, params);

# get Helmholtz operator
omega = 2 * pi * 20.0;
H = get_helmholtz_LU(params, omega);
b = zeros(Complex{Float64}, params.Nz * params.Nx);
idx = (src.isx+params.npml-1)*params.Nz + src.isz + params.ntop;
b[idx] = 100000.0

# snapshot
u = zeros(Complex{Float64}, params.Nz * params.Nx);
ldiv!(u, H, b);
u = reshape(u, params.Nz, params.Nx)

# plotting
SeisPlotTX(real(u), wbox=9, hbox=3);


# # generate recordings
# get_recordings!(rec, src, params, print_flag=true);
#
# # forward modeling in frequency domain
# pre = get_wavefield_FDFD(src, params; tmax=tmax, print_flag=true);

#plotting
# it = 400
# figure(figsize=(9,3)); imshow(pre_f[:,:,it*2], cmap="seismic");
# figure(figsize=(9,3)); imshow(pre_t[:,:,it], cmap="seismic");
#
# iz = 51; ix = 151;
# figure(); plot(pre_f[iz,ix,:]); plot(pre_t[iz,ix,:])
