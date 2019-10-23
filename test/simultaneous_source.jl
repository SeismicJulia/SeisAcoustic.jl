using SeisPlot, SeisAcoustic

dir_work = joinpath(homedir(), "Desktop/lsrtm");

path_rho = joinpath(dir_work, "physical_model/rho.rsf");
path_vel = joinpath(dir_work, "physical_model/vel.rsf");

# read the physical model
(hdr_rho, rho0) = read_RSdata(path_rho);
(hdr_vel, vel0) = read_RSdata(path_vel);

# # cropped model for imaging
vel = vel0[1:4:1350,650:4:2300];
rho = rho0[1:4:1350,650:4:2300];

SeisPlotTX(vel, hbox=3.38, wbox=4.13, cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
SeisPlotTX(rho, hbox=3.38, wbox=4.13, cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0007; tmax = 10.0;

# top boundary condition
free_surface = true;
data_format  = Float32;
order        = 5;

# tdparams for generating observations
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# initialize a source
isz = [2,2]; isx = [100, 300]; ot=[0.0, 5.0];
srcs= get_multi_sources(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");

# initialize recordings
irx = collect(1:2:fidiff.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, fidiff);

# save the wavefield
path_pre = joinpath(dir_work, "wavefield/pressure_100.rsf");
@time multi_step_forward!(rec, srcs, fidiff, path_pre=path_pre, interval=10);

# make animation
(hdr, d) = read_RSdata(path_wfd);
path_ani = joinpath(dir_work, "wavefield/animation_100");
SeisAnimation(path_ani, d; wbox=4.13, hbox=3.38, cmap="gray", vmin=-10, vmax=10, xticks=[], yticks=[], interval=20);
