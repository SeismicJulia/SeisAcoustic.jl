using SeisPlot, SeisAcoustic, PyPlot

dir_work = joinpath(homedir(), "Desktop");

path_rho = joinpath(dir_work, "physical_model/rho.rsf");
path_vel = joinpath(dir_work, "physical_model/vel.rsf");

# read the physical model
(hdr_rho, rho0) = read_RSdata(path_rho);
(hdr_vel, vel0) = read_RSdata(path_vel);

# # cropped model for imaging
vel = vel0[1:4:1350,150:4:2300];
rho = rho0[1:4:1350,150:4:2300];

# plotting velocity model
SeisPlotTX(vel, wbox=8, hbox=5, cmap="rainbow", vmax=maximum(vel), vmin=minimum(vel),
           ox=0, dx=0.00625, oy=0, dy=0.00625, xlabel="X (km)", ylabel="Z (km)");
tight_layout(); savefig("/Users/wenlei/Desktop/vel.pdf"); close();

# vertical and horizontal grid size
dz = 6.25; dx = 6.25;

# time step size and maximum modelling length
dt = 0.0007; tmax = 6.0;

# top boundary condition
free_surface = true;
data_format  = Float32;
order        = 5;

# tdparams for generating observations
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# receiver location
irx = collect(1: 2: fidiff.nx-10);
irz = 2 * ones(length(irx));

# source location at isx=100
isz = 2; isx = 90; ot=0.0;
src = Source(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");

# single source simulation
interval = 10;  # save one every 10 snapshots 
path_pre = joinpath(dir_work, "single_pre.rsf")
multi_step_forward!(src, fidiff, path_pre=path_pre, interval=interval);

# single source recording
rec = multi_step_forward!(irz, irx, src, fidiff);
a = quantile(abs.(rec.p[:]), (98/100));
SeisPlotTX(rec.p; wbox=4, hbox=5, cmap="seismic", vmax=a, vmin=-a,
           ox=1, dx=1, oy=0.0, dy=0.0007, xlabel="Traces", ylabel="Time (s)")
tight_layout(); savefig("/Users/wenlei/Desktop/rec.pdf"); close();


(hdr, p) = read_RSdata("/Users/wenlei/Desktop/single_pre.rsf");
a = quantile(abs.(p[:]), (98/100));
path_pre = "/Users/wenlei/Desktop/move_pre";
SeisAnimation(path_out, p; wbox=8, hbox=5, cmap="gray", vmax=a, vmin=-a,
              ox=0, dx=0.00625, oy=0, dy=0.00625, title="wavefield", xlabel="X (km)", ylabel="Z (km)", interval=70)

# make the record
p = rec.p[1:interval:end,:];
a = quantile(abs.(p[:]), (98/100));
(nt, nr) = size(p);
d = zeros(Float32, nt, nr, nt);
for i = 1 : nt
    d[1:i,:,i] .= p[1:i,:]
end
path_rec = "/Users/wenlei/Desktop/move_rec";
SeisAnimation(path_rec, d; wbox=4, hbox=5, cmap="gray", vmax=a, vmin=-a,
              ox=1, dx=1, oy=0.0, dy=0.007, xlabel="Traces", ylabel="Time (s)", interval=35)


# simultaneouse source survey
tmax = 12.0;
fidiff = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

isz = [2,2,2,2,2]; isx = [90, 180, 270, 360, 450]; ot=[0.0, 2.5, 3.1, 4.7, 5.3];
srcs= get_multi_sources(isz, isx, fidiff; ot=ot, amp=100000, fdom=20, type_flag="miniphase");


# single source simulation
interval = 20;
path_pre = joinpath(dir_work, "multiple_pre.rsf")
multi_step_forward!(srcs, fidiff, path_pre=path_pre, interval=interval);

# single source recording
rec = multi_step_forward!(irz, irx, srcs, fidiff);
a = quantile(abs.(rec.p[:]), (98/100));
SeisPlotTX(rec.p; wbox=4, hbox=10, cmap="gray", vmax=a, vmin=-a,
           ox=1, dx=1, oy=0.0, dy=0.0007, xlabel="Traces", ylabel="Time (s)")
tight_layout(); savefig("/Users/wenlei/Desktop/rec.pdf"); close();


(hdr, p) = read_RSdata("/Users/wenlei/Desktop/multiple_pre.rsf");
a = quantile(abs.(p[:]), (98/100));
path_pre = "/Users/wenlei/Desktop/multiple_pre";
SeisAnimation(path_pre, p; wbox=8, hbox=5, cmap="gray", vmax=a, vmin=-a,
              ox=0, dx=0.00625, oy=0, dy=0.00625, title="wavefield", xlabel="X (km)", ylabel="Z (km)", interval=70)

# make the record
p = rec.p[1:interval:end,:];
a = quantile(abs.(p[:]), (98/100));
(nt, nr) = size(p);
d = zeros(Float32, nt, nr, nt);
for i = 1 : nt
    d[1:i,:,i] .= p[1:i,:]
end
path_rec = "/Users/wenlei/Desktop/multiple_rec";
SeisAnimation(path_rec, d; wbox=4, hbox=10, cmap="gray", vmax=a, vmin=-a,
              ox=1, dx=1, oy=0.0, dy=0.014, xlabel="Traces", ylabel="Time (s)", interval=70)
