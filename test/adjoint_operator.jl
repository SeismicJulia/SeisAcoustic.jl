# homogeneous velocity and density model
nz = 101; nx = 201;
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 1.0;  # use second as unit

# receiver location
irx = collect(1 :2 :nx);
irz = 2 * ones(length(irx));

# ==============================================================================
#                   dot-product test for the adjoint wavefield
# ==============================================================================
test_type = [(true , Float32, 1.0e-6 ),
             (true , Float64, 1.0e-12),
             (false, Float32, 1.0e-6 ),
             (false, Float64, 1.0e-12)
            ]
@testset "multi-step forward modelling and its adjoint operator" begin

    for (free_surface, data_format, tol) in test_type

        # generate finite difference stencil
        params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                 data_format=data_format, fd_flag="taylor", order=3, npml=20, apml=900.);

        # allocate snapshots for forward and adjoint operator
        spt1_f = Snapshot(params);
        spt2_f = Snapshot(params);
        spt1_b = Snapshot(params);
        spt2_b = Snapshot(params);

        # initialize spt1_f and spt2_b with random numbers
        N = params.Nz * params.Nx
        spt1_f.vz .= randn(params.data_format, N); spt2_b.vz .= randn(params.data_format, N);
        spt1_f.vx .= randn(params.data_format, N); spt2_b.vx .= randn(params.data_format, N);
        spt1_f.pz .= randn(params.data_format, N); spt2_b.pz .= randn(params.data_format, N);
        spt1_f.px .= randn(params.data_format, N); spt2_b.px .= randn(params.data_format, N);

        # temporary variables
        tmp    = zeros(params.data_format, params.Nz * params.Nx);
        tmp_z1 = zeros(params.data_format, params.Nz);
        tmp_z2 = zeros(params.data_format, params.Nz);
        tmp_x1 = zeros(params.data_format, params.Nx);
        tmp_x2 = zeros(params.data_format, params.Nx);

        # nt-step forward
        nt = 78
        for it = 1 : nt
            one_step_forward!(spt2_f, spt1_f, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2);
            copy_snapshot!(spt1_f, spt2_f);
        end

        # nt-step adjoint
        for it = 1 : nt
            one_step_adjoint!(spt1_b, spt2_b, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2);
            copy_snapshot!(spt2_b, spt1_b);
        end

        # inner product
        tmp1 = (dot(spt1_f.vz, spt1_b.vz) + dot(spt1_f.vx, spt1_b.vx)
              + dot(spt1_f.pz, spt1_b.pz) + dot(spt1_f.px, spt1_b.px))

        tmp2 = (dot(spt2_f.vz, spt2_b.vz) + dot(spt2_f.vx, spt2_b.vx)
              + dot(spt2_f.pz, spt2_b.pz) + dot(spt2_f.px, spt2_b.px))

        @test (tmp1-tmp2) / tmp1 < tol
    end
end

# ==============================================================================
#           adjoint Pz = Px in computational area
# ==============================================================================
@testset "test pz == px for the adjoint wavefield" begin

    for (free_surface, data_format, tol) in test_type

        # generate finite difference stencil
        params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                 data_format=data_format, fd_flag="taylor", order=3, npml=20, apml=900.);


        # generate band-limited random recordings
        w   = ricker(12.0, params.dt)
        hl  = floor(Int64, length(w)/2)
        rec = Recordings(irz, irx, params);
        idx_o= [1]

        # generate data trace by trace
        for i = 1 : rec.nr
            tmp = conv(randn(params.nt)*1000, w)
            copyto!(rec.p, idx_o[1], tmp, hl+1, params.nt)
            idx_o[1] = idx_o[1] + params.nt
        end

        # save adjoint pz and px as an 3D cube
        path_spt = joinpath(homedir(), "Desktop/snapshot.rsf")
        multi_step_adjoint(rec, params; path_spt=path_spt)

        # # read the adjoint snapshot cube
        (hdr, d) = read_RSdata(path_spt);

        zl = params.ntop+1; zu = zl + params.nz - 1;
        xl = params.npml+1; xu = xl + params.nx - 1;

        @test norm(d[zl:zu,xl:xu,3,:] - d[zl:zu,xl:xu,4,:]) / norm(d[zl:zu,xl:xu,3,:]) < tol
    end
end
