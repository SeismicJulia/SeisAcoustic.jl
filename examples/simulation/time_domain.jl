using SeisPlot, SeisAcoustic

# velocity and density model
vel = 3000 * ones(100, 300);  # m/s
vel[51:end,:] .= 3500.0  # m/s
rho = 2000 * ones(100, 300);  # kg/m^3

# number of PML layers
npml = 20

# top boundary condition
free_surface = true    #(pml or free_surface)

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
# @time multi_step_forward!(rec, srcs, params);
# SeisPlotTX(rec.p, pclip=98)
#
# forward modeling and save the pressure field to hard drive
path_pre = joinpath(homedir(), "Desktop/pressure.rsb");
multi_step_forward!(path_pre, src, params);
(hdr, pre0) = read_RSdata(path_pre);

# save the boundary of wavefield, which can be used for reconstruction
path_bnd = joinpath(homedir(), "Desktop/boundary.rsb")
path_wfd = joinpath(homedir(), "Desktop/wavefield.rsb")
get_boundary_wavefield(path_bnd, path_wfd, src, params)

# read boundary and last wavefield
bnd  = read_boundary(path_bnd);
wfd  = read_one_wavefield(path_wfd, 1);

# reconstruct the pressure field forward
pre1 = pressure_reconstruct_forward(bnd, src, params);

# reconstruct the pressure field backward
pre2 = pressure_reconstruct_backward(bnd, wfd, src, params);


# path_pre = joinpath(homedir(), "Desktop/snashot.rsb");
# multi_step_forward!(path_pre, src, params; save_flag="snapshot");

path_spt = joinpath(homedir(), "Desktop/snashot.rsb");
multi_step_adjoint!(path_spt, rec, params; save_flag="snapshot")


pz = zeros(params.data_format, params.nz * params.nx);
px = zeros(params.data_format, params.nz * params.nx);
spt = read_one_snapshot(path_spt, 300);
for i = 1 : length(params.spt2wfd)
    j = params.spt2wfd[i]
    pz[i] = spt.pz[j]
    px[i] = spt.px[j]
end
pz = reshape(pz, params.nz, params.nx);
px = reshape(px, params.nz, params.nx);

figure(); imshow(pz, cmap="seismic")
figure(); imshow(px, cmap="seismic")

# # test the adjoint operator
# irx = collect(1:2:params.nx);
# irz = 2 * ones(length(irx));
# rec = Recordings(irz, irx, params);
# multi_step_forward!(rec, src, params)
#
# # the random inputs
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




# test the one-step adjoint operator
spt1_f = Snapshot(params);
spt2_f = Snapshot(params);

# initialize as random number
for ix = 1 : params.Nx
    row_idx = (ix-1) * params.Nz
    for iz = 1 : params.Nz
        idx= row_idx + iz
        spt1_f.vz[idx] = randn()
        spt1_f.vx[idx] = randn()
        spt1_f.pz[idx] = randn()
        spt1_f.px[idx] = randn()
    end
end

tmp_z1 = zeros(params.data_format, params.Nz);
tmp_z2 = zeros(params.data_format, params.Nz);
tmp_x1 = zeros(params.data_format, params.Nx);
tmp_x2 = zeros(params.data_format, params.Nx);

# one-step forward
one_step_forward!(spt2_f, spt1_f, params,
                  tmp_z1, tmp_z2, tmp_x1, tmp_x2);

# test the one-step adjoint operator
spt1_b = Snapshot(params);
spt2_b = Snapshot(params);

# initialize as random number
idx = 0
for ix = 1 : params.Nx
    row_idx = (ix-1) * params.Nz
    for iz = 1 : params.Nz
        idx= row_idx + iz
        spt1_b.vz[idx] = randn()
        spt1_b.vx[idx] = randn()
        spt1_b.pz[idx] = randn()
        spt1_b.px[idx] = spt1_b.pz[idx]
    end
end

tmp = randn(params.Nz * params.Nx);

# one-step forward
one_step_adjoint!(spt2_b, spt1_b, params,
                  tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

tmp0 = foo(spt1_f, spt2_b)
tmp2 = foo(spt2_f, spt1_b)

function foo(spt1, spt2)
    tmp = 0.
    tmp = tmp + dot(spt1.vz, spt2.vz)
    tmp = tmp + dot(spt1.vx, spt2.vx)
    tmp = tmp + dot(spt1.pz, spt2.pz)
    tmp = tmp + dot(spt1.px, spt2.px)
    return tmp
end

pz = zeros(params.nz * params.nx)
px = zeros(params.nz * params.nx)
for i = 1 : length(params.spt2wfd)
    j = params.spt2wfd[i]
    pz[i] = spt2_b.pz[j]
    px[i] = spt2_b.px[j]
end
pz = reshape(pz, params.nz, params.nx)
px = reshape(px, params.nz, params.nx)




















#
