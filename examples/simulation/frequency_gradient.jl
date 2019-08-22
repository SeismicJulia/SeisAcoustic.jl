using SeisPlot, SeisAcoustic, FFTW, LinearAlgebra

# velocity and density model
nz = 101; nx = 301
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = false    #(pml or free_surface)

# grid size
h = 10

# organize these parameters into a structure
fparams = FdParams(rho, vel, free_surface, h; data_format=Float64);

# generate observations of one frequency slice by a point source
isz = 5; isx = 151; dt = 0.001;
src = Source(isz, isx, dt, fparams; amp=1/rho[isz,isx]/vel[isz,isx]^2*10^5);

# total modelling length
tmax = 2.0;
nt   = floor(Int64, tmax/dt) + 1

# working frequency
df = 1.0 / (dt * nt);
f_work = 20.0;
iw = floor(Int64, f_work/df) + 1;
f_work = (iw-1) * df;
omega  = 2.0 * pi * f_work;
wlet   = zeros(nt)
copyto!(wlet, src.it_min, src.p)
fwave  = fft(wlet)

# LU decomposition of Helmhotz operator
H = get_helmholtz_LU(fparams, omega);

# source vector
s = zeros(Complex{Float64}, fparams.Nz*fparams.Nx)
s[src.src2spt] = fwave[iw]

# produce modelling true wave field
u = zeros(Complex{Float64}, fparams.Nz * fparams.Nx);
ldiv!(u, H, s);
SeisPlotTX(real(reshape(u, fparams.Nz, fparams.Nx)), wbox=9, hbox=3)

# specify the monochromatic recordings
irx = collect(1:2:nx);
irz = 2 * ones(Int64, length(irx));
dobs= MonochromaticRecordings(irz, irx, omega, fparams);
sample_spt2rec!(dobs, u);


# initial homogenous velocity model
vel = 3000 * ones(nz, nx)
fparams = FdParams(rho, vel, free_surface, h; data_format=Float64);

# helmholtz operator of homogeneous model
H = get_helmholtz_LU(fparams, omega);

# source side wavefield
u = zeros(Complex{Float64}, fparams.Nz * fparams.Nx);
ldiv!(u, H, s);
SeisPlotTX(real(reshape(u, fparams.Nz, fparams.Nx)), wbox=9, hbox=3)

dsyn = MonochromaticRecordings(irz, irx, omega, fparams);
sample_spt2rec!(dsyn, u);

# compute the residue
dres = get_residue(dsyn, dobs);
g = velocity_gradient(dres, u, H, fparams);
g = reshape(g, fparams.nz, fparams.nx)
SeisPlotTX(real(g), wbox=9, hbox=3)


# compute the gradient numerically
vel0    = 3000 * ones(nz, nx);
delta_m = 1e-7;

# one model parameter's gradient
iz = 37; ix = 137;

# positive velocity model partubation
vel0[iz,ix] = vel0[iz,ix] + delta_m

# model parameter
fparams = FdParams(rho, vel0, free_surface, h; data_format=Float64);
H = get_helmholtz_LU(fparams, omega);

# source side wavefield
u = zeros(Complex{Float64}, fparams.Nz * fparams.Nx);
ldiv!(u, H, s);

d1 = MonochromaticRecordings(irz, irx, omega, fparams);
sample_spt2rec!(d1, u);
SeisPlotTX(real(reshape(u, fparams.Nz, fparams.Nx)), wbox=9, hbox=3)

# compute the residue
dres = get_residue(d1, dobs);
J_plus = 1/2.0 * (norm(dres.p))^2

# negative velocity model partubation
vel0[iz,ix] = vel0[iz,ix] - 2*delta_m

fparams = FdParams(rho, vel0, free_surface, h; data_format=Float64);
H = get_helmholtz_LU(fparams, omega);

# source side wavefield
u = zeros(Complex{Float64}, fparams.Nz * fparams.Nx);
ldiv!(u, H, s);

d2 = MonochromaticRecordings(irz, irx, omega, fparams);
sample_spt2rec!(d2, u);

# compute the residue
dres = get_residue(d2, dobs);
J_minus = 1/2.0 * (norm(dres.p))^2;

# compute gradient numerically
g_num = (J_plus - J_minus) / (2 * delta_m)
g[iz,ix]
(g[iz,ix] - g_num) / g_num
