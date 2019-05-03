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

figure(figsize=(9,3)); imshow(pre0[:,:,200], cmap="seismic");
figure(figsize=(9,3)); imshow(pre1[:,:,200], cmap="seismic");
figure(figsize=(9,3)); imshow(pre2[:,:,200], cmap="seismic");

# true velocity model
vel1 = copy(vel);
vel1[51,151] = 3000*(1+1e-7);
params1 = ModelParams(rho, vel1, free_surface, dz, dx, dt, tmax;
          data_format=Float64, fd_flag="taylor", order=2, npml=20, apml=900.)

rec_obs = Recordings(irz, irx, params1);

# forward modeling of simultaneous sources
multi_step_forward!(rec_obs, src, params1);
SeisPlotTX(rec_obs.p, pclip=98);

# compute the residue
residue = Recordings(irz, irx, params);
for ir = 1 : rec_syn.nr
    for it = 1 : rec_syn.nt
        residue.p[it,ir] = rec_syn.p[it,ir] - rec_obs.p[it,ir]
    end
end
SeisPlotTX(residue.p, pclip=98);


# test the adjoint wavefield pz=px
path_adj = joinpath(homedir(), "Desktop/adjoint.rsb");
multi_step_adjoint!(path_adj, residue, params; save_flag="snapshot")

(hdr, d) = read_RSdata(path_adj);
d = reshape(d, hdr.n1*hdr.n2, hdr.n3, hdr.n4);
tmp= zeros(1)
for it = 1 : hdr.n4
    for i = 1 : length(params.spt2wfd)
        j = params.spt2wfd[i]
        tmp[1] = tmp[1] + (d[j,3,it] - d[j,4,it])^2
    end
    println("$it")
end

pz = zeros(params.nz * params.nx);
px = zeros(params.nz * params.nx);
it = 1
for i = 1 : length(params.spt2wfd)
    j = params.spt2wfd[i]
    pz[i] = d[j,3,it]
    px[i] = d[j,4,it]
end
pz = reshape(pz, params.nz, params.nx);
px = reshape(px, params.nz, params.nx);
SeisPlotTX(pz, wbox=6, hbox=2);
SeisPlotTX(px, wbox=6, hbox=2);




# compute the gradient with respect to velocity
@time g = bulk_gradient(residue, path_bnd, path_wfd, src, params);
g = reshape(g, params.nz, params.nx);
SeisPlotTX(g, pclip=98, wbox=9, hbox=3);
tmp = -g[:,150]; tmp[1:10] .= 0.0;
figure(); plot(tmp)
