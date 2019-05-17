using SeisPlot, SeisAcoustic

# ========================================================
#        testing the code to get source-side wave field
# ========================================================
# 2-layers velocity model
vel = 3000 * ones(101, 301);  # m/s
vel[51:end,:] .= 3500;

# constant density model
rho = 2000 * ones(101, 301);  # kg/m^3

# number of PML layers
npml = 20;

# top boundary condition
free_surface = true;

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);

# initialize a source
src = Source(2, 150, params; ot=0.0, fdom=20.0,
      type_flag="ricker", amp=100000, location_flag="index");

# initialize recordings
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec = Recordings(irz, irx, params);

# generate recordings and boundary value
path_bnd = joinpath(homedir(), "Desktop/bnd.rsb");
path_wfd = joinpath(homedir(), "Desktop/wfd.rsb");
multi_step_forward!(rec, src, params; path_bnd=path_bnd, path_wfd=path_wfd)
SeisPlotTX(rec.p)

# compute source-side wavefield
pre1 = source_side_wavefield(path_bnd, path_wfd, src, params);
pre2 = get_sourceside_wavefield(src, params; iflag=1);
pre3 = get_sourceside_wavefield(src, params; iflag=2);

# check the consistent of the value
norm(pre1 - pre2) / norm(pre1);
norm(pre1 - pre3) / norm(pre1);

it = 399
SeisPlotTX(pre1[:,:,it], cmap="seismic", wbox=9, hbox=3)
SeisPlotTX(pre2[:,:,it], cmap="seismic", wbox=9, hbox=3)
SeisPlotTX(pre3[:,:,it], cmap="seismic", wbox=9, hbox=3)

# ==============================================================================
#            check the analytical gradient against numerical gradient
# ==============================================================================
using SeisPlot, SeisAcoustic, PyPlot, LinearAlgebra

# size of model
nz = 101; nx = 301;

# 2-layers velocity model
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;

# constant density model
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = true;

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax; data_format=Float64);

# initialize a source
isz = 2; isx = 150;
src = Source(isz, isx, params; amp=100000);

# generate observed data
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
dobs= Recordings(irz, irx, params);

# generate recordings and boundary value
multi_step_forward!(dobs, src, params)
SeisPlotTX(dobs.p, cmap="seismic", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# ==============================================================================
#                      the gradient based on adjoint method
# initial velocity model
vel0 = 3000 * ones(nz, nx);
params0 = TdParams(rho, vel0, free_surface, dz, dx, dt, tmax; data_format=Float64);

# generate synthetic data and save the boundary value
path_bnd = joinpath(homedir(), "Desktop/bnd.rsb");
path_wfd = joinpath(homedir(), "Desktop/wfd.rsb");
dsyn = Recordings(irz, irx, params0);

multi_step_forward!(dsyn, src, params0, path_bnd=path_bnd, path_wfd=path_wfd);
SeisPlotTX(dsyn.p, cmap="seismic", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# compute the residue
dres = get_residue(dsyn, dobs);
SeisPlotTX(dres.p, cmap="seismic", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# get the gradient
@time g_adj = velocity_gradient(dres, path_bnd, path_wfd, src, params0);
g_adj = reshape(g_adj, nz, nx)
SeisPlotTX(g_adj, cmap="seismic", wbox=9, hbox=3)

# ==============================================================================
#                      compute gradient numerically
# initial velocity model as homogenous one
vel0    = 3000 * ones(nz, nx);
delta_m = 1e-7;

# one model parameter's gradient
iz = 51; ix = 201;

# positive velocity model partubation
vel0[iz,ix] = vel0[iz,ix] + delta_m

# model parameter
params0 = TdParams(rho, vel0, free_surface, dz, dx, dt, tmax; data_format=Float64);

# get synthetic data
dsyn = Recordings(irz, irx, params0);
multi_step_forward!(dsyn, src, params0);

# compute the residue
dres = get_residue(dsyn, dobs);
J_plus = 1/2.0 * (norm(dres.p))^2

# negative velocity model partubation
vel0[iz,ix] = vel0[iz,ix] - 2*delta_m

# model parameter
params0 = TdParams(rho, vel0, free_surface, dz, dx, dt, tmax; data_format=Float64);

# synthetic data
dsyn = Recordings(irz, irx, params0);
multi_step_forward!(dsyn, src, params0);

# compute the residue
dres = get_residue(dsyn, dobs);
J_minus = 1/2.0 * (norm(dres.p))^2;

# compute gradient numerically
g_num = (J_plus - J_minus) / (2 * delta_m)
g_adj[iz,ix]
g_adj[iz,ix] - g_num


# some codes used for testing
# multi_step_adjoint!(path_adj, dres, params0, save_flag="pressure")
# (hdr, adj) = read_RSdata(path_adj);
#
# # get the source side wave field
# pre = get_sourceside_wavefield(src, params0);
#
# # compute the gradient analytically
# g_ana = zeros(params.data_format, params.nz, params.nx);
# for it = 1 : params0.nt-1
#     for i2 = 1 : params.nx
#         for i1 = 1 : params.nz
#             g_ana[i1,i2] = g_ana[i1,i2] + pre[i1,i2,it] * adj[i1,i2,it+1] / 2.0
#         end
#     end
# end
#
# g_ana[iz,ix]
# g_num
# SeisPlotTX(g_ana, cmap="seismic", hbox=4, wbox=9)
