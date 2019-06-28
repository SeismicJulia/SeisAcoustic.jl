"""
   compute the gradient of objective function with respect to velocity model for one
frequency slice, r is the residue obtained by subtracting observed data from synthetic,
u is source side wavefield with PML boundary part, H is the LU factorization
of Helmhotz operator.
"""
function velocity_gradient_test(r::MonochromaticRecordings, u::Vector{Tv},
         H, params::FdParams) where {Tv<:Complex{Float64}}

    # insert the conplex conjugate of residue to snapshot
    b = Vector{Complex{Float64}}(undef, params.Nz*params.Nx)
    inject_rec2spt!(b, r)

    # square of radian frequency
    omega2 = r.omega * r.omega

    # some coefficients of mixed-grid
    w1 = 0.6287326
    w2 = 0.0928166
    w3 = (1.0 - w1 - 4.0*w2) / 4.0

    # create a table save the coefficients for mass acceleration average
    C  = OffsetArray{Float64}(undef, -1:1, -1:1)

    C[-1,-1] = w3; C[-1,0] = w2; C[-1,1] = w3;
    C[ 0,-1] = w2; C[ 0,0] = w1; C[ 0,1] = w2;
    C[ 1,-1] = w3; C[ 1,0] = w2; C[ 1,1] = w3;

    # compute adjoint wavefield
    u_adj = zeros(Complex{Float64}, params.Nz*params.Nx)
    ldiv!(u_adj, transpose(H), b)

    # apply the imaging condition to get the gradient
    g = zeros(Complex{Float64}, params.nz * params.nx)

    # both u and u_adj are vectors include PML layers
    for i2 = 1 : params.nx
        idx2 = (i2-1) * params.nz

        for i1 = 1 : params.nz
            idx1 = idx2 + i1

            # the constant part
            tmp   = params.spt2wfd[idx1]
            alpha = -2.0 * omega2 * u[tmp]
            alpha = alpha / (params.rho[i1,i2] * (params.vel[i1,i2])^3)

            for d2 = -1 : 1
                l2 = (i2+params.npml+d2-1) * params.Nz

                for d1 = -1 : 1
                    if 1 <= i1+d1 || !params.free_surface
                       l1 = l2 + i1 + d1 + params.ntop

                       g[idx1] = g[idx1] + C[d1,d2] * u_adj[l1]
                    end
                end
            end

            g[idx1] = g[idx1] * alpha
        end
    end

    return real(g)
end

using SeisPlot, PyPlot, SeisAcoustic, FFTW, OffsetArrays, LinearAlgebra

# velocity and density model
nz = 101; nx = 301
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# top boundary condition
free_surface = true    #(pml or free_surface)

# grid size
h = 10

# organize these parameters into a structure
fparams = FdParams(rho, vel, free_surface, h; data_format=Float64);

# generate observations of one frequency slice by a point source
isz = 5; isx = 151; dt = 0.001;
src = Source(isz, isx, dt, fparams);

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
g = velocity_gradient_test(dres, u, H, fparams);
SeisPlotTX(real(reshape(g, fparams.nz, fparams.nx)), wbox=9, hbox=3)


# compute the gradient numerically
vel0    = 3000 * ones(nz, nx);
delta_m = 1e-7;

# one model parameter's gradient
iz = 51; ix = 151;

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
g[iz,ix] - g_num

# some other test the transpose of LU decomposition
# A = sprand(Complex{Float64}, 10, 10, 0.3) + spdiagm(0=>rand(Complex{Float64}, 10));
# B = Matrix(transpose(A))
# B = Sparse(B
# H = lu(B)
# P = lu(A)
#
# x0 = rand(Complex{Float64}, 10);
# b  = Vector{Complex{Float64}}(undef, 10);
# mul!(b, transpose(A), x0);
#
# x = Vector{Complex{Float64}}(undef, 10);
# ldiv!(x, H, b)
#
# x1 = Vector{Complex{Float64}}(undef, 10);
# ldiv!(x1, transpose(P), b)
