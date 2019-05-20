"""
   immutable struct contain the model parameters for FDFD
*  data_format : Float32 or Float64
*  nz          : depth of model
*  nx          : length of model
*  npml        : number of PML layers
*  free_surface: true -> no obsorbtion on the surface
                 false-> PML on the surface
*  Nz          : depth of model (include PML padding)
*  Nx          : length of model (include PML padding)
*  ntop        : number of padded PML layers on the top
*  spt2wfd     : auxillary index vector for computing wavefield from snapshot
*  h           : grid spatial size
*  rho         : 2D density model
*  vel         : 2D P-wave velocity model
"""
struct FdParams{Ti<:Int64, Tv<:AbstractFloat}

    data_format  :: DataType
    nz           :: Ti
    nx           :: Ti
    npml         :: Ti
    apml         :: Tv
    free_surface :: Bool
    Nz           :: Ti
    Nx           :: Ti
    ntop         :: Ti
    spt2wfd      :: Vector{Ti}
    h            :: Tv
    rho          :: Matrix{Tv}
    vel          :: Matrix{Tv}
end

"""
   Overloading the show function for FDParams
"""
function show(io::IO, params::FdParams)

    # size of model
    @printf("nz = %4d, nx = %4d, npml = %4d\n", params.nz, params.nx, params.npml)
    @printf("h  = %4.1f\n", params.h)

    # top boundary condition
    if params.free_surface
       @printf("top boundary conditon = %20s\n", "free surface")
    else
       @printf("top boundary conditon = %20s\n", "PML obsorbing")
    end

    # data presicion (Float32, Float64)
    @printf("data_format=%16s\n", params.data_format)

    return nothing
end

"""
   get the auxillary index vectors which facilitate extracting wavefield from snapshot
"""
function get_spt2wfd_index(nz::Ti, nx::Ti, npml::Ti, free_surface::Bool) where {Ti<:Int64}

    # size of wavefield
    len = nz * nx

    # allocate 3 vectors of integer
    index = zeros(Int64, len)

    # size of padded model
    if free_surface
       Nz = nz + npml
       ntop = 0
    else
       Nz = nz + npml * 2
       ntop = npml
    end
    Nx = nx + npml * 2

    # map snapshot to wavefield
    icount = 0
    for i2 = 1 : nx
        tmp= (i2 + npml - 1)*Nz + ntop

        for i1 = 1 : nz
            icount = icount + 1
            index[icount] = tmp + i1
        end
    end

    return index
end

"""
   create the parameter structure for frequency domain finite difference method
"""
function FdParams(rho, vel, free_surface::Bool, h;
         data_format=Float32, npml=20, apml=900.)

    # get the size of the model
    (nz, nx) = size(rho)
    if (nz, nx) != size(vel)
       error("model size mismatch")
    end

    # unify the data format
    rho = convert(Matrix{data_format}, rho)
    vel = convert(Matrix{data_format}, vel)
    h   = convert(data_format, h)
    apml= convert(data_format, apml)

    # size of model after padding PML layers
    if free_surface
       Nz   = nz +     npml
       ntop = 0
    else
       Nz   = nz + 2 * npml
       ntop = npml
    end
    Nx = nx + 2 * npml

    # compute auxillary index vectors
    idx_spt2wfd = get_spt2wfd_index(nz, nx, npml, free_surface)

    # call the default struct constructor
    return FdParams(data_format, nz, nx, npml, apml, free_surface, Nz, Nx, ntop, idx_spt2wfd, h, rho, vel)
end

"""
   initialize the pml damping profile in z- and x-direction
"""
function DampingProfile(Nz::Ti, Nx::Ti, npml::Ti, apml::Tv,
         free_surface::Bool, h::Tv, omega::Float64) where {Ti<:Int64, Tv<:AbstractFloat}

    # convert the data format of radian frequency
    apml  = convert(Float64, apml )
    h     = convert(Float64, h    )

    # get the damping profile in Z direction
    if free_surface
       (zdamp_integer, zdamp_stagger) = get_fsdamping_profile(Nz, npml, apml, h, omega)
    else
       (zdamp_integer, zdamp_stagger) = get_damping_profile(Nz, npml, apml, h, omega)
    end

    # get the damping profile in X direction
    (xdamp_integer, xdamp_stagger) = get_damping_profile(Nx, npml, apml, h, omega)

    # # data format convert
    # zdamp_integer = convert(Vector{Complex{params.data_format}}, zdamp_integer)
    # zdamp_stagger = convert(Vector{Complex{params.data_format}}, zdamp_stagger)
    # xdamp_integer = convert(Vector{Complex{params.data_format}}, xdamp_integer)
    # xdamp_stagger = convert(Vector{Complex{params.data_format}}, xdamp_stagger)

    # group into a struct
    return zdamp_integer, zdamp_stagger, xdamp_integer, xdamp_stagger
end

"""
   set up damping profile for PML boundary conditions, n is the model size include pml,
npml is the number of PML layers, omega is radian frequency.
"""
function get_damping_profile(n::Ti, npml::Ti, apml::Tv, h::Tv, omega::Tv) where {Ti<:Int64, Tv<:Float64}

    # constants
    pi_over_2 = 0.5 * pi   # pi/2

    # damping profile of PML
    damp_integer = ones(Complex{Float64}, n+2)
    damp_stagger = ones(Complex{Float64}, n+2)

    # PML width
    pml_width = npml * h

    # model width include PML layers
    model_width = (n - 1) * h

    # both side of damp_integer and one side of damp_stagger
    for i = 1 : npml

        # current coordinate of x
        x_integer = (i-1) * h
        x_stagger = (i-1+0.5) * h

        # distance of current grid to the edge of computation area
        dis_integer = pml_width - x_integer
        dis_stagger = pml_width - x_stagger

        # pml coefficients
        gamma_integer = apml * (1.0 - cos(dis_integer / pml_width * pi_over_2))
        gamma_stagger = apml * (1.0 - cos(dis_stagger / pml_width * pi_over_2))

        # dampping profile (one extra point) 1/(1-iγ/w)
        damp_integer[i+1] = 1.0 / (1.0 - im * gamma_integer / omega)
        damp_stagger[i+1] = 1.0 / (1.0 - im * gamma_stagger / omega)

        # symmetrical damp_integer
        damp_integer[n-i+2] = damp_integer[i+1]
    end

    # one extra point facilitate building Helmhotz operator
    damp_integer[1]   = damp_integer[2]
    damp_integer[n+2] = damp_integer[n+1]

    # right side of damp_stagger
    for i = 1 : npml + 1

        # the coordinate of the stagger grid (from right to left)
        x_stagger = model_width + 0.5 * h - (i-1)*h

        # distance of current grid to the edge of computation area
        dis_stagger = x_stagger - (model_width - pml_width)

        # pml coefficients
        gamma_stagger = apml * (1.0 - cos(dis_stagger / pml_width * pi_over_2))

        # damping profile
        damp_stagger[n-i+2] = 1.0 / (1.0 - im * gamma_stagger / omega)
    end

    damp_stagger[1]   = damp_stagger[2]
    damp_stagger[n+2] = damp_stagger[n+1]

    return damp_integer, damp_stagger
end

"""
   get the damping profile for vertical direction when free surface boundary
condition is implemented for the top boundary.
"""
function get_fsdamping_profile(Nz::Ti, npml::Ti, apml::Tv, h::Tv, omega::Tv) where {Ti<:Int64, Tv<:Float64}

    # constants
    pi_over_2 = 0.5 * pi   # pi/2

    # damping profile of PML
    damp_integer = ones(Complex{Float64}, Nz+2)
    damp_stagger = ones(Complex{Float64}, Nz+2)

    # PML width
    pml_width = npml * h

    # model width include PML layers
    model_width = (Nz - 1) * h

    # the bottom side
    for i = 1 : npml + 1

        # the coordinate of the stagger grid (from right to left)
        x_integer = model_width           - (i-1)*h
        x_stagger = model_width + 0.5 * h - (i-1)*h

        # distance of current grid to the edge of computation area
        dis_integer = x_integer - (model_width - pml_width)
        dis_stagger = x_stagger - (model_width - pml_width)

        # pml coefficients
        gamma_integer = apml * (1.0 - cos(dis_integer / pml_width * pi_over_2))
        gamma_stagger = apml * (1.0 - cos(dis_stagger / pml_width * pi_over_2))

        # damping profile
        damp_integer[Nz-i+2] = 1.0 / (1.0 - im * gamma_integer / omega)
        damp_stagger[Nz-i+2] = 1.0 / (1.0 - im * gamma_stagger / omega)
    end

    damp_integer[Nz+2] = damp_integer[Nz+1]
    damp_stagger[Nz+2] = damp_stagger[Nz+1]

    return damp_integer, damp_stagger
end

"""
   compute complex bulk modulus which support visco-acoustic wave-equation, padding
one extra layer. The bulk modulus is padded with one more outer layer.
"""
function vel2bulk(rho::Matrix{Tv}, vel::Matrix{Tv}) where {Tv<:Float64}

    # size of model
    (nz, nx)    = size(rho)

    # allocate space for bulk modulus
    bulk = zeros(Float64, nz+2, nx+2)
    for ix = 1 : nx
        for iz = 1 : nz
            bulk[iz+1,ix+1] = rho[iz,ix] * vel[iz,ix] * vel[iz,ix]
        end
    end

    # visco-acoustic case
    # bulk = zeros(Complex{data_format}, nz+2, nx+2)
    # for ix = 1 : nx
    #     for iz = 1 : nz
    #         vc = vel[iz,ix] * (1.0 - 0.5*im / q[iz,ix])
    #         mu[iz+1,ix+1] = rho[iz,ix] * vc * vc
    #     end
    # end

    # left and right
    for iz = 1 : nz
        bulk[iz+1, 1]    = bulk[iz+1, 2]
        bulk[iz+1, nx+2] = bulk[iz+1, nx+1]
    end

    # upper and bottom
    for ix = 1 : nx
        bulk[1, ix+1]    = bulk[2, ix+1]
        bulk[nz+2, ix+1] = bulk[nz+1,ix+1]
    end

    # four corners
    bulk[1, 1]       = bulk[2, 2]
    bulk[1, nx+2]    = bulk[2, nx+1]
    bulk[nz+2, 1   ] = bulk[nz+1, 2]
    bulk[nz+2, nx+2] = bulk[nz+1, nx+1]

    return bulk
end

"""
   compute buoyancy from density, the computed buoyancy is paded with one extra layer.
"""
function rho2buoyancy(rho::Matrix{Tv}) where {Tv <: Float64}

    # size of the model
    (nz, nx) = size(rho)

    # allocate space for buoyancy
    buoyancy = zeros(Float64, nz+2, nx+2)
    for ix = 1 : nx
        for iz = 1 : nz
            if rho[iz, ix] != 0.0
               buoyancy[iz+1, ix+1] = 1.0 / rho[iz,ix]
            else
               error("density is zero")
            end
        end
    end

    # left and right layer
    for iz = 1 : nz
        buoyancy[iz+1, 1]    = buoyancy[iz+1, 2]
        buoyancy[iz+1, nx+2] = buoyancy[iz+1, nx+1]
    end

    # top and bottom layer
    for ix = 1 : nx
        buoyancy[1, ix+1]    = buoyancy[2, ix+1]
        buoyancy[nz+2, ix+1] = buoyancy[nz+1, ix+1]
    end

    # four corners
    buoyancy[1, 1]       = buoyancy[2, 2   ]
    buoyancy[1, nx+2]    = buoyancy[2, nx+1]
    buoyancy[nz+2, 1]    = buoyancy[nz+1, 2]
    buoyancy[nz+2, nx+2] = buoyancy[nz+1, nx+1]

    return buoyancy
end

"""
   get the buoyancy at stagger grid by averaging the buoyancy at reference grid
"""
function buoyancy_average!(bu::OffsetArray, buoyancy::Matrix{Tv}, i::Ti, j::Ti) where {Ti<:Int64, Tv<:Float64}

    # some constants
    one_over_two  =  1.0 / 2.0
    one_over_four =  1.0 / 4.0

    # first column
    bu[-1,-1] = one_over_four * (buoyancy[i-1,j-1] + buoyancy[i,j-1] + buoyancy[i-1,j] + buoyancy[i,j])
    bu[0 ,-1] = one_over_two  * (buoyancy[i,j-1]   + buoyancy[i,j]                                    )
    bu[1 ,-1] = one_over_four * (buoyancy[i,j-1] + buoyancy[i+1,j-1] + buoyancy[i,j] + buoyancy[i+1,j])

    # second column
    bu[-1,0] = one_over_two * (buoyancy[i,j] + buoyancy[i-1,j])
    bu[0 ,0] =                 buoyancy[i,j]
    bu[1 ,0] = one_over_two * (buoyancy[i,j] + buoyancy[i+1,j])

    # third column
    bu[-1,1] = one_over_four * (buoyancy[i-1,j] + buoyancy[i,j] + buoyancy[i-1,j+1] + buoyancy[i,j+1])
    bu[0 ,1] = one_over_two  * (buoyancy[i,j] + buoyancy[i,j+1]                                      )
    bu[1 ,1] = one_over_four * (buoyancy[i,j] + buoyancy[i+1,j] + buoyancy[i,j+1] + buoyancy[i+1,j+1])

    return nothing
end

"""
   the kernal of building Helmhotz operator, output the LU factorization of Helmhotz operator
"""
function get_helmholtz_LU(params::FdParams, omega)

    # convert radian frequency to Float64
    omega = convert(Float64, omega)

    # constants for mass acceleration averaging
    w1 = 0.6287326  # c
    w2 = 0.0928166  # d
    w3 = (1.0 - w1 - 4.0*w2) / 4.0 # (1-c-4d)/4

    # weights for stagger and rotated grid
    alpha = 0.5617365     # α
    beta  = 1.0 - alpha   # 1-α

    # omega square and one over h square
    omega2       = omega * omega
    one_over_h2  = 1.0  / (params.h * params.h)
    one_over_4h2 = 0.25 / (params.h * params.h)

    # model parameter padding to accomdate PML layers
    vel = model_padding(params.vel, params.npml, params.free_surface)
    rho = model_padding(params.rho, params.npml, params.free_surface)

    # convert the model parameter to Float64
    if params.data_format == Float32
       rho = convert(Matrix{Float64}, rho)
       vel = convert(Matrix{Float64}, vel)
    end

    # add support for visco-acoustic in the future
    # qfactor = model_padding(params.qfactor, params.npml, params.free_surface)

    # compute bulk modulus
    k        = vel2bulk(rho, vel)
    buoyancy = rho2buoyancy(rho)

    # compute the damping profile
    (zdamp_integer, zdamp_stagger, xdamp_integer, xdamp_stagger) = DampingProfile(params.Nz, params.Nx,
                                                                   params.npml, params.apml, params.free_surface, params.h, omega)

    # count the number of non-zero elements
    nnz = 0
    if params.free_surface
       nnz = nnz + 1 + 1 + 4 * 2              # four corners
       nnz = nnz + params.Nx - 2              # upper
       nnz = nnz + (params.Nx-2) * 6          # bottom
    else
       nnz = nnz + 4*4                        # four corners
       nnz = nnz + (params.Nx-2)*6*2          # upper and bottom
    end
    nnz = nnz + (params.Nz-2)*6*2             # left and right
    nnz = nnz + (params.Nz-2)*(params.Nx-2)*9 # inner part

    # allocate space for sparse Helmhotz matrix in coordinate format
    row_idx = zeros(Int64, nnz) # each row has 9 non-zeros elements
    col_idx = zeros(Int64, nnz)
    nzval   = zeros(Complex{Float64}, nnz)
    icount  = 0

    # buoyancy at stagger grid
    b = OffsetArray{Float64}(undef, -1:1, -1:1)

    for i2 = 1 : params.Nx
        for i1 = 1 : params.Nz

            # the corresponding index for buoyancy and bulk modulus as one extra more layer is padded
            i = i1 + 1
            j = i2 + 1

            # buoyancy at stagger grid
            buoyancy_average!(b, buoyancy, i, j)

            # damping parameter in z direction
            d1o = zdamp_integer[i]
            d1p = zdamp_stagger[i]
            d1m = zdamp_stagger[i-1]

            # damping parameter in x direction
            d2o = xdamp_integer[j]
            d2p = xdamp_stagger[j]
            d2m = xdamp_stagger[j-1]

            # row index
            row_val = (i2-1) * params.Nz + i1

            # node (0, 0)
            col_val = (i2-1  ) * params.Nz + i1
            icount  = icount + 1
            row_idx[icount] = row_val
            col_idx[icount] = col_val
            if (i1 == 1) && params.free_surface
               nzval[icount] = omega2 / k[i,j]
            else
               nzval[icount] = (-w1 * omega2 / k[i,j]
                                +alpha * one_over_h2  * ( d1o * (b[1,0]*d1p + b[-1,0]*d1m)
                                                         +d2o * (b[0,1]*d2p + b[0,-1]*d2m))
                                +beta  * one_over_4h2 * ( d1o * (b[1,1]*d1p + b[-1,-1]*d1m + b[1,-1]*d1p + b[-1,1]*d1m)
                                                         +d2o * (b[1,1]*d2p + b[-1,-1]*d2m + b[-1,1]*d2p + b[1,-1]*d2m)))
            end

            # node(-1,-1)
            if (i1 > 1) && (i2 > 1)
               col_val = (i2-1-1) * params.Nz + i1-1
               icount  = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount] = (-w3   * omega2 / k[i-1,j-1]
                                -beta * one_over_4h2 * (b[-1,-1]*d1o*d1m + b[-1,-1]*d2o*d2m))
            end

            # node(0,-1)
            if (i2 > 1)
               if (i1 != 1) || !params.free_surface
                  col_val = (i2-1-1) * params.Nz + i1
                  icount  = icount + 1
                  row_idx[icount] = row_val
                  col_idx[icount] = col_val
                  nzval[icount]   = (-w2 * omega2 / k[i,j-1]
                                     -alpha  * one_over_h2 * b[0,-1] * d2o * d2m
                                     -beta   * one_over_4h2* ( d2o*(b[-1,-1]*d2m + b[1,-1]*d2m)
                                                              -d1o*(b[-1,-1]*d1m + b[1,-1]*d1p)))
               end
            end

            # node(1,-1)
            if (i1 < params.Nz) && (i2 > 1)
               if (i1 != 1) || !params.free_surface
                  col_val = (i2-1-1) * params.Nz + i1+1
                  icount  = icount + 1
                  row_idx[icount] = row_val
                  col_idx[icount] = col_val
                  nzval[icount]   = (-w3   * omega2 / k[i+1,j-1]
                                     -beta * one_over_4h2 * (b[1,-1]*d1o*d1p + b[1,-1]*d2o*d2m))
               end
            end

            # node(-1, 0)
            if (i1 > 1)
               col_val = (i2-1  ) * params.Nz + i1-1
               icount = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount]   = (-w2 * omega2 / k[i-1,j]
                                  -alpha * one_over_h2  * b[-1,0] * d1o * d1m
                                  -beta  * one_over_4h2 * ( d1o*(b[-1,-1]*d1m + b[-1,1]*d1m)
                                                           -d2o*(b[-1,-1]*d2m + b[-1,1]*d2p)))
            end

            # node(1, 0)
            if (i1 < params.Nz)
               if (i1 != 1) || !params.free_surface
                  col_val = (i2-1  ) * params.Nz + i1+1
                  icount  = icount + 1
                  row_idx[icount] = row_val
                  col_idx[icount] = col_val
                  nzval[icount]   = (-w2 * omega2 / k[i+1,j]
                                     -alpha * one_over_h2  * b[1,0] * d1o * d1p
                                     -beta  * one_over_4h2 * ( d1o*(b[1,1]*d1p + b[1,-1]*d1p)
                                                              -d2o*(b[1,1]*d2p + b[1,-1]*d2m)))
               end
            end

            # node(-1, 1)
            if (i1 > 1) && (i2 < params.Nx)
               col_val = (i2-1+1) * params.Nz + i1-1
               icount  = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount]   = (-w3 * omega2 / k[i-1,j+1]
                                  -beta * one_over_4h2 * (b[-1,1]*d1o*d1m + b[-1,1]*d2o*d2p))
            end

            # node (0, 1)
            if (i2 < params.Nx)
               if (i1 != 1) || !params.free_surface
                  col_val = (i2-1+1) * params.Nz + i1
                  icount  = icount + 1
                  row_idx[icount] = row_val
                  col_idx[icount] = col_val
                  nzval[icount]   = (-w2 * omega2 / k[i,j+1]
                                     -alpha * one_over_h2  * b[0,1]*d2o*d2p
                                     -beta  * one_over_4h2 * ( d2o*(b[-1,1]*d2p + b[1,1]*d2p)
                                                              -d1o*(b[1,1]*d1p + b[-1,1]*d1m)))
               end

            end

            # node(1, 1)
            if (i1 < params.Nz) && (i2 < params.Nx)
               if (i1 != 1) || !params.free_surface
                  col_val = (i2-1+1) * params.Nz + i1+1
                  icount  = icount + 1
                  row_idx[icount] = row_val
                  col_idx[icount] = col_val
                  nzval[icount]   = (-w3 * omega2 / k[i+1,j+1]
                                     -beta * one_over_4h2 * (b[1,1]*d1o*d1p + b[1,1]*d2o*d2p))
               end
            end

        end # i1
    end # i2

    # check the size
    if nnz != icount
       error("number of non-zero element is wrong")
    end

    # save the sparse matrix in CSC format
    H = sparse(row_idx, col_idx, nzval)

    # return the LU factorization of Helmhotz operator
    return lu(H)
end

"""
   generate source for frequency domain finite difference
"""
function Source(sz, sx, params::FdParams; location_flag="index", ot=0.0, dt=0.001, amp=1.0, fdom=20.0, hl=128,
                type_flag="ricker", p=Vector{Float32}(undef,0))

    # source location given index
    if location_flag == "index"
       isz = round(Int64, sz)
       isx = round(Int64, sx)

    # source location given as distance
    elseif location_flag == "distance"
       isz = round(Int64, sz/params.h) + 1
       isx = round(Int64, sx/params.h) + 1
    else
       error("wrong specification of source location")
    end

    # error checking
    if isz > params.nz || isx > params.nx || isz < 1 || isx < 1
       error("source located outside of modeling area")
    end

    # can't put on the free-surface
    if params.free_surface && isz == 1
       error("can't inject source on the surface")
    end

    # the auxillary index mapping the location of source to snapshot or wavefield
    src2spt = (isx+params.npml-1) * params.Nz + isz + params.ntop
    src2wfd = (isx            -1) * params.nz + isz

    # non-user provided source wavelet
    if length(p) == 0
       if type_flag == "ricker"
          p = ricker(fdom, dt)
          p = amp * p
       elseif type_flag == "miniphase"
          p = ricker(fdom, dt)
          p = convert2miniphase(p)
          p = amp * p
       elseif type_flag == "sinc"
          fc = params.data_format(fdom)
          hl = round(Int64, hl)
          p  = tapered_sinc(fc, hl, dt)
          p  = amp * p
       end
    end

    # make the data type consistent
    p = convert(Vector{params.data_format}, p)
    nt= length(p)

    # the index time range of source wavelet
    it_min = floor(Int64, ot/dt) + 1
    it_max = it_min + nt - 1

    # call the default constructor
    return Source(isz, isx, src2spt, src2wfd, dt, it_min, it_max, p)
end

"""
   get the forward modeling wavefield via frequency domain finite difference
"""
function get_wavefield_FDFD(src::Source, params::FdParams;
         flower=0.0, fupper=60.0, tmax=2.0, print_flag=false)

    # length of wavefield
    N = params.nz * params.nx

    # number of time samples
    nt = floor(Int64, tmax / src.dt) + 1

    # allocate space for 3D wavefield [nz*nx, nt]
    pressure = zeros(Complex{Float64}, N, nt)

    # frequency interval
    fmax = 1.0  / src.dt
    df   = fmax / nt

    # lower frequency index bound
    iw_lower = floor(Int64, flower / df) + 1
    if iw_lower < 2
       iw_lower = 2
    end

    # upper frequency index bound
    iw_upper = ceil(Int64, fupper / df) + 1
    if iw_upper > floor(Int64, nt / 2) + 1
       iw_upper = floor(Int64, nt / 2) + 1
       fmax = (iw_upper-1)*df
    end

    # allocate space for source time function
    wlet = zeros(nt)
    copyto!(wlet, src.it_min, convert(Vector{Float64},src.p))

    # transform to frequency domain
    fwave = fft(wlet)

    # allocate space for left and right hand side
    b = zeros(Complex{Float64}, params.Nz * params.Nx)
    u = zeros(Complex{Float64}, params.Nz * params.Nx)

    # loop over all frequency slice
    for iw = iw_lower : iw_upper

        # radian frequency
        omega = 2.0 * pi * (iw-1) * df

        # build damping profile
        H = get_helmholtz_LU(params, omega)

        #the correct way to implement source
        bulk = params.rho[src.isz, src.isx] * params.vel[src.isz, src.isx] * params.vel[src.isz, src.isx]
        b[src.src2spt] = im * omega * fwave[iw] / bulk
        ldiv!(u, H, b)

        # save the pressure component
        for i = 1 : N
            tmp = u[params.spt2wfd[i]]
            pressure[i, iw] = tmp
            pressure[i, nt-iw+2] = conj(tmp)
        end

        # print the frequency slice
        if print_flag
           println("finished $iw")
        end
    end

    # transform back to time domain
    pressure = real(ifft(pressure, 2))

    # return a 3D cube
    return reshape(pressure, params.nz, params.nx, nt)
end

"""
   constructor for recordings for frequency domain finite difference
"""
function Recordings(rz::Vector, rx::Vector, dt, tmax,
         params::FdParams; location_flag="index")

    # number of receivers
    nr  = length(rz)
    if length(rx) != nr
       error("length of receiver coordinate vector doesn't match")
    end

    irz = zeros(Int64, nr)
    irx = zeros(Int64, nr)

    # location are provided as index
    if location_flag == "index"
       for i = 1 : nr
           irz[i] = round(Int64, rz[i])
           irx[i] = round(Int64, rx[i])
       end

    # location are given as distance
    elseif location_flag == "distance"
       for i = 1 : nr
           irz[i] = round(Int64, rz[i]/params.h) + 1
           irx[i] = round(Int64, rx[i]/params.h) + 1
       end
    else
       error("wrong specification of receiver location")
    end

    # error checking
    for i = 1 : nr
        if irz[i] > params.nz || irx[i] > params.nx || irz[i] < 1 || irx[i] < 1
           error("receiver located outside of modeling area")
        end
    end

    # can't put receivers on the free surface
    if params.free_surface
       for i = 1 : nr
           if irz[i] == 1
              error("can't put receiver at free surface, move it deeper")
           end
       end
    end

    # convert dt format
    dt = convert(params.data_format, dt)
    nt = floor(Int64, tmax/dt) + 1

    # the auxillary vector mapping snapshot to recordings
    spt2rec = zeros(Int64, nr)
    for i = 1 : nr
        spt2rec[i] = (irx[i]+params.npml-1) * params.Nz + irz[i] + params.ntop
    end

    # initalize empty recordings
    rec= Recordings(nt, nr, dt, irz, irx, spt2rec,
         zeros(params.data_format, nt, nr))

    return rec
end

"""
   get the recordings via frequency domain finite difference
"""
function get_recordings!(rec::Recordings, src::Source, params::FdParams;
         flower=0.0, fupper=60.0, print_flag=false) where {Ti <: Int64}

    # length of wavefield
    N = params.nz * params.nx

    # number of time samples
    nt = rec.nt

    # check the time sample interval
    if rec.dt != src.dt
       error("time sample interval doesn't match")
    end

    # frequency interval
    fmax = 1.0  / src.dt
    df   = fmax / nt

    # lower frequency index bound
    iw_lower = floor(Int64, flower / df) + 1
    if iw_lower < 2
       iw_lower = 2
    end

    # upper frequency index bound
    iw_upper = ceil(Int64, fupper / df) + 1
    if iw_upper > floor(Int64, nt / 2) + 1
       iw_upper = floor(Int64, nt / 2) + 1
       fmax = (iw_upper-1)*df
    end

    # allocate space for source time function
    wlet = zeros(nt)
    copyto!(wlet, src.it_min, convert(Vector{Float64},src.p))

    # transform to frequency domain
    fwave = fft(wlet)

    # allocate space for left and right hand side
    b = zeros(Complex{Float64}, params.Nz * params.Nx)
    u = zeros(Complex{Float64}, params.Nz * params.Nx)

    # the pressure component
    p = zeros(Complex{Float64}, nt, rec.nr)

    # loop over all frequency slice
    for iw = iw_lower : iw_upper

        # radian frequency
        omega = 2.0 * pi * (iw-1) * df

        # build damping profile
        H = get_helmholtz_LU(params, omega)

        #the correct way to implement source
        bulk = params.rho[src.isz, src.isx] * params.vel[src.isz, src.isx] * params.vel[src.isz, src.isx]
        b[src.src2spt] = im * omega * fwave[iw] / bulk
        ldiv!(u, H, b)

        # save the pressure component
        for i = 1 : rec.nr
            idx = rec.spt2rec[i]
            p[iw,i]     = u[idx]
            p[nt-iw+2,i]= conj(u[idx])
        end

        # print the frequency slice
        if print_flag
           println("finished $iw")
        end
    end

    # transform back to time domain
    ifft!(p, 1)
    p = real(p)

    # copy to recordings
    copyto!(rec.p, 1, p)

    return nothing
end

# vel = model_padding(params.vel, params.npml, params.free_surface)
# rho = model_padding(params.rho, params.npml, params.free_surface)
# k        = vel2bulk(rho, vel)
# buoyancy = rho2buoyancy(rho)

# """
#    convert sparse matrix from coordinate format to CSC format
# """
# function COOCSC(m::Ti, n::Ti,
#          irow::Vector{Ti}, icol::Vector{Ti}, nzval::Vector{Tv}) where {Ti<:Int64, Tv<:Number}
#
#   etype = eltype(nzval)
#   nnz   = length(nzval)
#
#   colptr = zeros(Int64, n+1)
#   rowidx = zeros(Int64, nnz)
#   nzval1 = zeros(etype, nnz)
#
#   # Determine the column length
#   for i = 1 : nnz
#       colptr[icol[i]] = colptr[icol[i]] + 1
#   end
#
#   # starting position of the element in each column
#   k = 1
#   for j = 1 : n + 1
#       k0 = colptr[j]
#       colptr[j] = k
#       k = k + k0
#   end
#
#   # fill in the output
#   for k = 1 : nnz
#       i = irow[k]
#       j = icol[k]
#       x = nzval[k]
#       idx = colptr[j]
#       nzval1[idx] = x
#       rowidx[idx] = i
#       colptr[j  ] = idx + 1
#   end
#
#   # shift back colptr
#   for j = n : -1 : 1
#       colptr[j+1] = colptr[j]
#   end
#   colptr[1] = 1
#
#   A = spzeros(m, n)
#   A.rowval = rowidx
#   A.colptr = colptr
#   A.nzval  = nzval1
#
#   return A
#
# end
#
# """
#    For sparse matrix in coordinate form, remove the repeated ones
# """
# function remove_repeated_samples(m::Ti, n::Ti,
#          irow::Ti, icol::Ti, nzval::Tv) where {Ti<:Int64, Tv<:Number}
#
#   etype = eltype(nzval)
#   nnz   = length(nzval)
#
#   colptr = zeros(Int64, n+1)
#   rowidx = zeros(Int64, nnz)
#   nzval1 = zeros(etype, nnz)
#
#   # Determine the column length
#   for i = 1 : nnz
#       colptr[icol[i]] = colptr[icol[i]] + 1
#   end
#
#   # starting position of each column
#   k = 1
#   for j = 1 : N + 1
#       k0 = colptr[j]
#       colptr[j] = k
#       k = k + k0
#   end
#
#   # fill in the output
#   for k = 1 : nnz
#       i = irow[k]
#       j = icol[k]
#       x = nzval[k]
#       idx = colptr[j]
#       nzval1[idx] = x
#       rowidx[idx] = i
#       colptr[j  ] = idx + 1
#   end
#
#   # shift back colptr
#   for j = N : -1 : 1
#       colptr[j+1] = colptr[j]
#   end
#   colptr[1] = 1
#
#   A = spzeros(m, n)
#   A.rowval = rowidx
#   A.colptr = colptr
#   A.nzval  = nzval1
#
#   return A
#
# end
