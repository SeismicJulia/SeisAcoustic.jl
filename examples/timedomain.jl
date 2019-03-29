using SeisAcoustic

vel = 3000 * ones(100, 300);
rho = ones(100, 300);

npml = 20
free_surface = false

dz = 10; dx = 10;

dt = 0.001; tmax = 2.0;
fdom = 10;

params = ModelParams(rho, vel, npml, free_surface, dz, dx, dt, tmax, fdom; data_format=Float32)
ofds = ObsorbFDStencil(params);
rfds
