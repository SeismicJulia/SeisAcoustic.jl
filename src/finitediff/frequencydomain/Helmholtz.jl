"""
   damping profile for 2D frequency domain finite difference
"""
struct DampingProfile{Tv<:AbstractFloat}
    zdamp_integer :: Vector{Complex{Tv}}
    zdamp_stagger :: Vector{Complex{Tv}}
    xdamp_integer :: Vector{Complex{Tv}}
    xdamp_stagger :: Vector{Complex{Tv}}
end

"""
   initialize the pml damping profile in z- and x-direction
"""
function DampingProfile(params::ModelParams, omega::Tv;
         apml=900.0) where {Tv<:AbstractFloat}

    # this code only handle equal grid size
    if params.dz != params.dx
       error("the grid size doesn't equal")
    end

    # get the damping profile in Z direction
    (zdamp_integer, zdamp_stagger) = get_damping_profile(params.Nz, params.npml, params.dz, omega; apml=apml)

    # get the damping profile in X direction
    (xdamp_integer, xdamp_stagger) = get_damping_profile(params.Nx, params.npml, params.dx, omega; apml=apml)

    # data format convert
    zdamp_integer = convert(Vector{Complex{params.data_format}}, zdamp_integer)
    zdamp_stagger = convert(Vector{Complex{params.data_format}}, zdamp_stagger)
    xdamp_integer = convert(Vector{Complex{params.data_format}}, xdamp_integer)
    xdamp_stagger = convert(Vector{Complex{params.data_format}}, xdamp_stagger)

    # group into a struct
    return DampingProfile(zdamp_integer, zdamp_stagger, xdamp_integer, xdamp_stagger)
end

"""
   set up damping profile for PML boundary conditions, n is the model size include pml,
npml is the number of PML layers, omega is radian frequency.
"""
function get_damping_profile(n::Ti, npml::Ti, h, omega; apml=900.0) where {Ti<:Int64}

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
   compute complex bulk modulus which support visco-acoustic wave-equation, padding
one extra layer.
"""
function vel2bulk(rho::Matrix{Tv}, vel::Matrix{Tv}) where {Tv<:AbstractFloat}

    # data format
    data_format = eltype(rho)

    # size of model
    (nz, nx)    = size(rho)

    # allocate space for bulk modulus
    bulk = zeros(data_format, nz+2, nx+2)
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
   compute buoyancy from density, padding one extra layer.
"""
function rho2buoyancy(rho::Matrix{Tv}) where {Tv<: AbstractFloat}

    data_format = eltype(rho)
    (nz, nx)    = size(rho)

    # allocate space for buoyancy
    buoyancy = zeros(data_format, nz+2, nx+2)
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
function buoyancy_average!(bu::OffsetArray, buoyancy::Matrix{Tv}, i::Ti, j::Ti) where {Ti<:Int64, Tv<:AbstractFloat}

    data_format = eltype(buoyancy)

    # some constants
    one_over_two  = data_format(1.0 / 2.0)
    one_over_four = data_format(1.0 / 4.0)

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
function helmholtz_kernal(k::Matrix{Tv}, buoyancy::Matrix{Tv}, damp::DampingProfile,
         omega::Tv, params::ModelParams) where {Tv<:AbstractFloat}

    # constants for mass acceleration averaging
    w1 = params.data_format(0.6287326)  # c
    w2 = params.data_format(0.0928166)  # d
    w3 = params.data_format((1.0 - w1 - 4.0*w2) / 4.0) # (1-c-4d)/4

    # weights for stagger and rotated grid
    alpha = params.data_format(0.5617365  )    # α
    beta  = params.data_format(1.0 - alpha)    # 1-α

    # omega square and one over h square
    omega2       = params.data_format(omega * omega)
    one_over_h2  = params.data_format(1.0  / (params.dz * params.dz))
    one_over_4h2 = params.data_format(0.25 / (params.dz * params.dz))

    # number of non-zero elements
    nnz = 0
    nnz = nnz + 4*4                           # four corners
    nnz = nnz + (params.Nx-2)*6*2             # upper and bottom
    nnz = nnz + (params.Nz-2)*6*2             # left and right
    nnz = nnz + (params.Nz-2)*(params.Nx-2)*9 # inner part

    # allocate space for sparse Helmhotz matrix in coordinate format
    row_idx = zeros(Int64, nnz) # each row has 9 non-zeros elements
    col_idx = zeros(Int64, nnz)
    nzval   = zeros(Complex{params.data_format}, nnz)
    icount  = 0

    # buoyancy at stagger grid
    b = OffsetArray{params.data_format}(undef, -1:1, -1:1)

    for i2 = 1 : params.Nx
        for i1 = 1 : params.Nz

            # the corresponding index for buoyancy and bulk modulus
            i = i1 + 1
            j = i2 + 1

            # buoyancy at stagger grid
            buoyancy_average!(b, buoyancy, i, i)

            # damping parameter in z direction
            d1o = damp.zdamp_integer[i]
            d1p = damp.zdamp_stagger[i]
            d1m = damp.zdamp_stagger[i-1]

            # damping parameter in x direction
            d2o = damp.xdamp_integer[j]
            d2p = damp.xdamp_stagger[j]
            d2m = damp.xdamp_stagger[j-1]

            # row index
            row_val = (i2-1) * params.Nz + i1

            # node (0, 0)
            col_val = (i2-1  ) * params.Nz + i1
            icount  = icount + 1
            row_idx[icount] = row_val
            col_idx[icount] = col_val
            nzval[icount] = (-w1 * omega2 / k[i,j]
                             +alpha * one_over_h2  * ( d1o * (b[1,0]*d1p + b[-1,0]*d1m)
                                                      +d2o * (b[0,1]*d2p + b[0,-1]*d2m))
                             +beta  * one_over_4h2 * ( d1o * (b[1,1]*d1p + b[-1,-1]*d1m + b[1,-1]*d1p + b[-1,1]*d1m)
                                                      +d2o * (b[1,1]*d2p + b[-1,-1]*d2m + b[-1,1]*d2p + b[1,-1]*d2m)))

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
               col_val = (i2-1-1) * params.Nz + i1
               icount  = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount]   = (-w2 * omega2 / k[i,j-1]
                                  -alpha  * one_over_h2 * b[0,-1] * d2o * d2m
                                  -beta   * one_over_4h2* ( d2o*(b[-1,-1]*d2m + b[1,-1]*d2m)
                                                           -d1o*(b[-1,-1]*d1m + b[1,-1]*d1p)))
            end

            # node(1,-1)
            if (i1 < params.Nz) && (i2 > 1)
               col_val = (i2-1-1) * params.Nz + i1+1
               icount  = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount]   = (-w3   * omega2 / k[i+1,j-1]
                                  -beta * one_over_4h2 * (b[1,-1]*d1o*d1p + b[1,-1]*d2o*d2m))
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
               col_val = (i2-1  ) * params.Nz + i1+1
               icount  = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount]   = (-w2 * omega2 / k[i+1,j]
                                  -alpha * one_over_h2  * b[1,0] * d1o * d1p
                                  -beta  * one_over_4h2 * ( d1o*(b[1,1]*d1p + b[1,-1]*d1p)
                                                           -d2o*(b[1,1]*d2p + b[1,-1]*d2m)))
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
               col_val = (i2-1+1) * params.Nz + i1
               icount  = icount + 1
               row_idx[icount] = row_val
               col_idx[icount] = col_val
               nzval[icount]   = (-w2 * omega2 / k[i,j+1]
                                  -alpha * one_over_h2  * b[0,1]*d2o*d2p
                                  -beta  * one_over_4h2 * ( d2o*(b[-1,1]*d2p + b[1,1]*d2p)
                                                           -d1o*(b[1,1]*d1p + b[-1,1]*d1m)))
             end

             # node(1, 1)
             if (i1 < params.Nz) && (i2 < params.Nx)
                 col_val = (i2-1+1) * params.Nz + i1+1
                 icount  = icount + 1
                 row_idx[icount] = row_val
                 col_idx[icount] = col_val
                 nzval[icount]   = (-w3 * omega2 / k[i+1,j+1]
                                    -beta * one_over_4h2 * (b[1,1]*d1o*d1p + b[1,1]*d2o*d2p))
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
   get the Helmholtz operator for one radian frequency and return its LU decomposition.
"""
function get_helmholtz_LU(params::ModelParams, omega)

    # make data format consistent
    omega = params.data_format(omega)

    # only support obsorbing surface boundary condition
    if params.free_surface == true
       error("only support obsorbing surface boundary condition")
    end

    # model parameter padding to accomdate PML layers
    vel = model_padding(params.vel, params.npml, params.free_surface)
    rho = model_padding(params.rho, params.npml, params.free_surface)

    # add support for visco-acoustic in future
    # qfactor = model_padding(params.qfactor, params.npml, params.free_surface)

    # compute bulk modulus
    k        = vel2bulk(rho, vel)
    buoyancy = rho2buoyancy(rho)

    # compute the damping profile
    damp = DampingProfile(params, omega)

    # compute helmhotz operator
    return helmholtz_kernal(k, buoyancy, damp, omega, params)
end

"""
   get the forward modeling wavefield via frequency domain finite difference
"""
function get_wavefield_FDFD(src::Source, params::ModelParams;
         flower=0.0, fupper=60.0, print_flag=false)

    # allocate space for 3D wavefield (nz, nx, nt)
    N = params.nz * params.nx
    pressure = zeros(Complex{Float64}, N, params.nt)

    # frequency interval
    fmax = 1.0  / params.dt
    df   = fmax / params.nt

    # lower frequency index bound
    iw_lower = floor(Int64, flower / df) + 1
    if iw_lower < 2
       iw_lower = 2
       @warn "can't compute when frequency <= 0Hz, minimum frequency starts from $df"
    end

    # upper frequency index bound
    iw_upper = ceil(Int64, fupper / df) + 1
    if iw_upper > floor(Int64, params.nt/2) + 1
       iw_upper = floor(Int64, params.nt/2) + 1
       fmax = (iw_upper-1)*df
       @warn "the maximum frequency must less than nyquist frequency, maximum frequency ends at $fmax"
    end

    # allocate space for source time function
    wlet = zeros(params.nt)
    copyto!(wlet, src.it_min, src.p)

    # transform to frequency domain
    fwave = fft(wlet) / sqrt(params.nt)

    # allocate space for left and right hand side
    b = zeros(Complex{Float64}, params.Nz * params.Nx)
    u = zeros(Complex{Float64}, params.Nz * params.Nx)

    # only support obsorbing surface boundary condition
    if params.free_surface == true
       error("only support obsorbing surface boundary condition")
    end

    # padding model parameter to accomdate PML layers
    vel = model_padding(params.vel, params.npml, params.free_surface)
    rho = model_padding(params.rho, params.npml, params.free_surface)
    k        = vel2bulk(rho, vel)
    buoyancy = rho2buoyancy(rho)

    # loop over all frequency slice
    for iw = iw_lower : iw_upper

        # radian frequency
        omega = 2.0 * pi * (iw-1) * df

        # build damping profile
        damp = DampingProfile(params, omega)
        H = helmholtz_kernal(k, buoyancy, damp, omega, params)

        #the correct way to implement source
        b[src.src2spt] = im * omega * fwave[iw] / k[src.src2spt]
        ldiv!(u, H, b)

        # save the pressure component
        for i = 1 : N
            tmp = u[params.spt2wfd[i]]
            pressure[i, iw] = tmp
            pressure[i, params.nt-iw+2] = conj(tmp)
        end

        # print the frequency slice
        if print_flag
           println("finished $iw")
        end
    end

    # transform back to time domain
    pressure = real(ifft(pressure, 2)) * sqrt(params.nt)

    # return a 3D cube
    return reshape(pressure, params.nz, params.nx, params.nt)
end

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
