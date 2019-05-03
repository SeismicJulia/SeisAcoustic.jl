using SeisPlot, SeisAcoustic

# homogeneous velocity and density model
vel = 3000 * ones(101, 301);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(101, 301);  # kg/m^3

# number of PML layers
npml = 20;

# top boundary condition
free_surface = true;   #(pml or free_surface)

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;  # use second as unit

# organize these parameters into a structure
params = ModelParams(rho, vel, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);
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
# srcs = get_multi_sources(isz, isx, params; amp=100000, ot=ot, fdom=15);

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
pre = source_side_wavefield(path_bnd, path_wfd, src, params);
SeisPlotTX(pre[:,:,200], cmap="seismic", wbox=9, hbox=3);

path_sfd = joinpath(homedir(), "Desktop/sfd.rsb");
hdr = RegularSampleHeader(pre, title="source-side wavefield");
write_RSdata(path_sfd, hdr, pre);

path_sfd1 = joinpath(homedir(), "Desktop/sfd_forward.rsb");
get_sourceside_wavefield(path_sfd1, src, params);

path_sfd2 = joinpath(homedir(), "Desktop/sfd_ground.rsb");
get_wavefield(path_sfd2, src, params);

(hdr, pre0) = read_RSdata(path_sfd);
(hdr, pre1) = read_RSdata(path_sfd1);
(hdr, pre2) = read_RSdata(path_sfd2);
norm(pre0 - pre1) / norm(pre0);
norm(pre0 - pre2) / norm(pre0);

it = 399
SeisPlotTX(pre0[:,:,it], cmap="seismic", wbox=9, hbox=3)
SeisPlotTX(pre1[:,:,it], cmap="seismic", wbox=9, hbox=3)
SeisPlotTX(pre2[:,:,it], cmap="seismic", wbox=9, hbox=3)


# ==============================================================================
#                       check the correctness of gradient
# ==============================================================================
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
rec_obs = Recordings(irz, irx, params);

# generate observations
multi_step_forward!(rec_obs, src, params);
SeisPlotTX(rec_obs.p, pclip=98, cmap="seismic");

# initial velocity model as a homogeneous one
vel1    = 3000 * ones(101, 301);  # m/s
delta_m = 1e-7;            # velocity pertubation
iz = 51; ix = 151;
vel1[iz, ix] = vel1[iz,ix] + delta_m

params1 = ModelParams(rho, vel1, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);

rec_syn = Recordings(irz, irx, params1);
multi_step_forward!(rec_syn, src, params1);

# compute the residue
residue = Recordings(irz, irx, params1);
for ir = 1 : rec_syn.nr
    for it = 1 : rec_syn.nt
        residue.p[it,ir] = rec_syn.p[it,ir] - rec_obs.p[it,ir]
    end
end
J_plus = norm(residue.p) / 2


# pertubate velocity
vel1[iz, ix] = vel1[iz,ix] - 2*delta_m

params1 = ModelParams(rho, vel1, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);

rec_syn = Recordings(irz, irx, params1);
multi_step_forward!(rec_syn, src, params1);

# compute the residue
residue = Recordings(irz, irx, params1);
for ir = 1 : rec_syn.nr
    for it = 1 : rec_syn.nt
        residue.p[it,ir] = rec_syn.p[it,ir] - rec_obs.p[it,ir]
    end
end
J_minus = norm(residue.p) / 2

# compute gradient numerically
g_num = (J_plus - J_minus) / (2 * delta_m)

# compute the gradient analytically
vel1    = 3000 * ones(101, 301);  # m/s

params1 = ModelParams(rho, vel1, free_surface, dz, dx, dt, tmax;
         data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.);

rec_syn = Recordings(irz, irx, params1);
multi_step_forward!(rec_syn, src, params1);

# compute the residue
residue = Recordings(irz, irx, params1);
for ir = 1 : rec_syn.nr
    for it = 1 : rec_syn.nt
        residue.p[it,ir] = rec_syn.p[it,ir] - rec_obs.p[it,ir]
    end
end

# save the adjoint wavefield
path_adj = joinpath(homedir(), "Desktop/adj.rsb");
multi_step_adjoint!(path_adj, residue, params1, save_flag="snapshot")
(hdr, adj) = read_RSdata(path_adj)

path_sfd = joinpath(homedir(), "Desktop/sfd_ground.rsb");
get_wavefield(path_sfd, src, params1)
(hdr, sfd) = read_RSdata(path_sfd)

path_vz = joinpath(homedir(), "Desktop/dvdz.rsb");
path_vx = joinpath(homedir(), "Desktop/dvdx.rsb");
get_dvdzx(path_vz, path_vx, src, params1)
(hdr, svz) = read_RSdata(path_vz)
(hdr, svx) = read_RSdata(path_vx)

g_ana = zeros(1)
i1 = iz + params1.ntop
i2 = ix + params1.npml
for it = 1 : params.nt-1
    g_ana[1] = g_ana[1] + svz[i1, i2, it] * adj[i1, i2, 3, it+1] + svx[i1, i2, it] * adj[i1, i2, 4, it+1]
end
g_ana[1] = g_ana[1] * rho[iz,ix] * vel1[iz,ix] * 2
