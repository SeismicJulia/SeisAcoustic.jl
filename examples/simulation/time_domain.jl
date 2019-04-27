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
params = ModelParams(rho, vel, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.)
# shows the default value for the keyword parameters
# data_format = (Float32 or Float64)
# fd_flag     = ("taylor" or "ls")
# order       = 2 - 10 if we use "ls" to compute the FD coefficients
#             = 2 - n  if we use "taylor" expansion to compute FD coefficients

# initialize a source
src = Source(2, 150, params; ot=0.0, fdom=20.0,
      type_flag="ricker", amp=100000, location_flag="index");

# initialize multi-sources
# isx = collect(5:60:295); ns=length(isx); isz = 2*ones(ns);
# ot  = 0.5*rand(ns);
# srcs = get_multi_sources(isz, isx, params; amp=100000, ot=ot);

# initialize recordings
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, params);

# forward modeling of simultaneous sources
@time multi_step_forward!(rec, src, params);
# @time multi_step_forward!(rec, src , ofds, params);
# imshow(rec.p, cmap="seismic", aspect=0.1)

# forward modeling and save the pressure field to hard drive
path_pre = joinpath(homedir(), "Desktop/pressure.rsb");
multi_step_forward(path_pre, src, ofds, params);
(hdr, pre0) = read_RSdata(path_pre);

# # save the boundary of wavefield, which can be used for reconstruction
path_bnd = joinpath(homedir(), "Desktop/boundary.rsb")
path_wfd = joinpath(homedir(), "Desktop/wavefield.rsb")
get_boundary_wavefield(path_bnd, path_wfd, src, ofds, params)

# read boundary and last wavefield
bnd  = read_boundary(path_bnd);
wfd  = read_one_wavefield(path_wfd, 1);

# reconstruct the pressure field forward
pre1 = pressure_reconstruct_forward(bnd, rfds, src, params);

# reconstruct the pressure field backward
pre2 = pressure_reconstruct_backward(bnd, wfd, rfds, src, params);


# test the adjoint operator
# irx = collect(1:2:params.nx);
# irz = 1 * ones(length(irx));
# rec = Recordings(irz, irx, params);
# multi_step_forward!(rec, src, ofds, params)
#
#
# # the random input
# w    = ricker(20.0, params.dt)
# hl   = floor(Int64, length(w)/2)
#
# rec1 = Recordings(irz, irx, params);
# idx_o= [1]
# for i = 1 : rec1.nr
#     tmp = conv(randn(params.nt)*1000, w)
#     copyto!(rec1.p, idx_o[1], tmp, hl+1, params.nt)
#     idx_o[1] = idx_o[1] + params.nt
# end
# @time p = multi_step_adjoint(rec1, ofds, src, params);
#
# tmp1 = dot(p, src.p)
# tmp2 = dot(vec(rec.p), vec(rec1.p))








#
