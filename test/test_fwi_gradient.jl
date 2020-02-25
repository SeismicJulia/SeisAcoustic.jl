using SeisPlot, SeisAcoustic, LinearAlgebra

# define true velocity model
nz = 101; nx = 301;
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;
rho = 2000 * ones(nz, nx);  # kg/m^3

# plot the physical model
SeisPlotTX(vel, wbox=9, hbox=3);

# number of PML layers
npml = 20;

# free surface for top boundary condition
free_surface = true;
dz = 10   ; dx   = 10 ;
dt = 0.001; tmax = 2.0;

params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);


# source term
src = Source(2, 150, params; ot=0.0, fdom=20.0,
      type_flag="ricker", amp=100000, location_flag="index")

# receiver
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, params);

# true recordings
multi_step_forward!(rec, src, params);
