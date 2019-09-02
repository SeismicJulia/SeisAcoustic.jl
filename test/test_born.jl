# ==============================================================================
#            check the analytical gradient against numerical gradient
# ==============================================================================
using SeisPlot, SeisAcoustic, LinearAlgebra, DSP

# size of model
nz = 101; nx = 301;

# 2-layers velocity model
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;

# constant density model
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = false;

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 2.0;

# precision
data_format=Float64;
order=3

# organize these parameters into a structure
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                  data_format=data_format, order=order);

# initialize a source
isz = 2; isx = 150;
src = Source(isz, isx, params; amp=100000, fdom=15);

# isx = collect(5:10:295); ns=length(isx);
# isz = 2   * ones(ns);
# ot  = 0.5 * rand(ns);
# pol = rand(Bool, ns);
# amp = 100000 * ones(ns);
# for i = 1 : ns
#     if !pol[i]
#        amp[i] = - amp[i]
#     end
# end
# srcs = get_multi_sources(isz, isx, params; amp=amp, ot=ot, fdom=5);

# generate observed data
irx = collect(1:2:params.nx);
irz = 2 * ones(length(irx));
dobs= Recordings(irz, irx, params);

# generate recordings and boundary value
@time multi_step_forward!(dobs, src, params);
SeisPlotTX(dobs.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# ==============================================================================
#                      the gradient based on adjoint method
# initial velocity model
vel0 = 3000 * ones(nz, nx);
params0 = TdParams(rho, vel0, free_surface, dz, dx, dt, tmax;
                   data_format=data_format, order=order);

# generate synthetic data and save the boundary value
path_bnd = joinpath(homedir(), "Desktop/bnd.rsb");
path_wfd = joinpath(homedir(), "Desktop/wfd.rsb");
dsyn = Recordings(irz, irx, params0);

# the source-side wavefield is generated simultanes source
isx = collect(5:10:295); ns=length(isx);
isz = 2   * ones(ns);
ot  = 0.5 * rand(ns);
pol = rand(Bool, ns);
amp = 100000 * ones(ns);
for i = 1 : ns
    if !pol[i]
       amp[i] = - amp[i]
    end
end
srcs = get_multi_sources(isz, isx, params; amp=amp, ot=ot, fdom=5);

# born approximation forward
multi_step_forward!(dsyn, srcs, params0, path_bnd=path_bnd, path_wfd=path_wfd);
SeisPlotTX(dsyn.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# compute the residue
dres = get_residue(dsyn, dobs);
SeisPlotTX(dres.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# set the scatter
# m = zeros(data_format, params.nz, params.nx)
# m[50,:] .= 100000.;
# m = vec(m);

# born approximation forward
m1 = 10000. * randn(params.nz * params.nx);
rec1 = Recordings(irz, irx, params0);
@time born_approximation_forward!(rec1, m1, path_bnd, srcs, params0);
SeisPlotTX(rec1.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# make band-limited random recordings
w    = ricker(15.0, params.dt)
hl   = floor(Int64, length(w)/2)
rec2 = Recordings(irz, irx, params);
idx_o= [1]
for i = 1 : rec2.nr
    tmp = conv(randn(params0.nt)*1000, w)
    copyto!(rec2.p, idx_o[1], tmp, hl+1, params0.nt)
    idx_o[1] = idx_o[1] + params0.nt
end
SeisPlotTX(rec2.p, cmap="gray", pclip=98,
           xlabel="Traces", ylabel="Time (s)", dy=0.001)

# born approximation adjoint
m2 = born_approximation_adjoint(rec2, path_bnd, path_wfd, srcs, params0);

# dot product test
tmp1 = dot(m1, m2)

p1 = rec1.p[2:end,:];
p2 = rec2.p[2:end,:];
tmp2 = dot(vec(p1), vec(p2))

(tmp1-tmp2) / tmp1


# examine the result
# i = 50
# p1 = dres.p[:,i] / maximum(abs.(dres.p[:,i]))
# p2 = d_born.p[:,i] / maximum(abs.(d_born.p[:,i]))
# figure(); plot(p1); plot(p2);
