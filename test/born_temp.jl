using SeisAcoustic

# homogeneous velocity and density model
nz = 51; nx = 75;
vel = 3000 * ones(nz, nx);  # m/s
vel[26:end,:] .= 3500;      # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 0.2;  # use second as unit

# receiver location
irx = collect(1 :2 :nx);
irz = 2 * ones(Int64, length(irx));

# source location
isx1 = floor(Int64, nx/2); isz1 = 2;
isx2 = collect(5:15:nx-5); ns=length(isx2); isz2 = 2*ones(ns);

test_type = [(true , Float32, 1.0e-6 ),
             (true , Float64, 1.0e-12),
             (false, Float32, 1.0e-6 ),
             (false, Float64, 1.0e-12)
            ]

(free_surface, data_format, tol) = test_type[1]


# generate finite difference stencil
params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
         data_format=data_format, fd_flag="taylor", order=3, npml=20, apml=900.);

# create multiple source term
source_term = get_multi_sources(isz2, isx2, params; ot=0.0*rand(ns), fdom=15.0, amp=100000);

# pre-compute source side wavefield and their boundary
dir_sourceside = joinpath(homedir(), "Desktop/source_side");
get_wavefield_bound(dir_sourceside, source_term, params);

# forward born approximation
# generate input for forward operator
m1      = randn(params.data_format, params.nz*params.nx)
path_m1 = joinpath(homedir(), "Desktop/m1.rsf");
hdr1    = RegularSampleHeader(m1);
write_RSdata(path_m1, hdr1, m1);

dir_born1 = joinpath(homedir(), "Desktop/forward_born");
born_approximation_forward(dir_born1, path_m1, irz, irx, dir_sourceside, params);


# adjoint operator
# adjoint born approximation
w   = ricker(15.0, params.dt)
hl  = floor(Int64, length(w)/2)

dir_born2 = joinpath(homedir(), "Desktop/adjoint_born");
rm(dir_born2, force=true, recursive=true)
mkdir(dir_born2)
for i = 1 : length(source_term)
    # generate data trace by trace
    rec2 = Recordings(irz, irx, params);
    for i = 1 : rec2.nr
        tmp = conv(randn(params.nt)*1000, w)
        tmp = convert(Vector{params.data_format}, tmp[hl+1:end-hl])
        rec2.p[:,i] .= tmp
    end

    file_name = join(["recordings_" "$i" ".bin"])
    path_rec  = joinpath(dir_born2, file_name)

    write_recordings(path_rec, rec2)
end
