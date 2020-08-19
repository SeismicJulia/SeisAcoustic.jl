# homogeneous velocity and density model
nz = 101; nx = 201;
vel = 3000 * ones(nz, nx);  # m/s
vel[51:end,:] .= 3500;  # m/s
rho = 2000 * ones(nz, nx);  # kg/m^3

# vertical and horizontal grid size
dz = 10; dx = 10;

# time step size and maximum modelling length
dt = 0.001; tmax = 1.0;  # use second as unit

# source location
isx1 = floor(Int64, nx/2); isz1 = 2;
isx2 = collect(5:40:nx-5); ns=length(isx2); isz2 = 2*ones(ns);

test_type = [(isz1, isx1, true , Float32, 1.0e-5 ),
             (isz1, isx1, true , Float64, 1.0e-10),
             (isz1, isx1, false, Float32, 1.0e-5 ),
             (isz1, isx1, false, Float64, 1.0e-10),
             (isz2, isx2, true , Float32, 1.0e-5 ),
             (isz2, isx2, true , Float64, 1.0e-10),
             (isz2, isx2, false, Float32, 1.0e-5 ),
             (isz2, isx2, false, Float64, 1.0e-10),
            ]
@testset "forward modelling and pressure field reconstruction (forward and backward)" begin

    for (isz, isx, free_surface, data_format, tol) in test_type

        # generate finite difference stencil
        params = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                 data_format=data_format, fd_flag="taylor", order=3, npml=20, apml=900.);

        # create source term
        ns = length(isz)
        if ns == 1
           source_term = Source(isz, isx, params; ot=0.0, fdom=15.0, amp=100000);

        elseif ns > 1
           source_term = get_multi_sources(isz, isx, params; ot=0.3*rand(ns), fdom=15.0, amp=100000);
        end

        # save the pressure field and boundary value
        path_pre  = joinpath(homedir(), "Desktop/pressure.rsf");
        path_bnd  = joinpath(homedir(), "Desktop/boundary_value.rsf");
        path_lwfd = joinpath(homedir(), "Desktop/last_wfd.rsf");

        # forward modeling of simultaneous sources
        multi_step_forward(source_term, params; path_pre=path_pre, path_bnd=path_bnd, path_lwfd=path_lwfd);
        (hdr, pre1) = read_RSdata(path_pre);

        # reconstruct pressure field using boundary infomation
        pre2 = pressure_reconstruct_forward(path_bnd, source_term, params);

        # reconstruct pressure field using boundary infomation
        pre3 = pressure_reconstruct_backward(path_bnd, path_lwfd, source_term, params);

        # compute the difference
        @test sqrt(sum((pre1 .- pre2) .^ 2)) / sqrt(sum(pre1 .^ 2)) < tol
        @test sqrt(sum((pre1 .- pre3) .^ 2)) / sqrt(sum(pre1 .^ 2)) < tol
    end
end
