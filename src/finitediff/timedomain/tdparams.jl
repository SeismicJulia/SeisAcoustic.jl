"""
   immutable struct contain the model parameters
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
*  spt2bnd     : auxillary index vector for computing wavefield boundaries from snapshot
*  bnd2wfd     : auxillary index vector for inserting wavefield boundaries to wavefield
*  dz          : vertical grid size
*  dx          : horizontal grid size
*  dt          : time step size
*  tmax        : maximum modelling length
*  nt          : number of time steps
*  rho         : 2D density model
*  vel         : 2D P-wave velocity model
*  MvzBvz      : FD stencil for updating v_z from v_z, vector length of Nz
*  MvzBp       : FD stencil for updating v_z from p  , vector length of Nz * Nx
*  dpdz        : sparse matrix computing dp/dz

*  MvxBvx      : FD stencil for updating v_x from v_x, vector length of Nx
*  MvxBp       : FD stencil for updating v_x from p  , vector length of Nz * Nx
*  dpdx        : sparse matrix computing dp/dx

*  MpzBpz      : FD stencil for updating p_z from p_z, vector length of Nz
*  MpzBvz      : FD stencil for updating p_z from v_z, vector length of Nz * Nx
*  dvdz        : sparse matrix computing dvz/dz

*  MpxBpx      : FD stencil for updating p_x from p_x, vector length of Nx
*  MpxBvx      : FD stencil for updating p_x from v_x, vector length of Nz * Nx
*  dvdx        : sparse matrix computing dvx/dx

*  RvzBp       : FD stencil for updating v_z from p, vector length of nz * nx
*  rpdz        : sparse matrix computing dp/dz

*  RvxBp       : FD stencil for updating v_x from p, vector length of nz * nx
*  rpdx        : sparse matrix computing dp/dx

*  RpzBvz      : FD stencil for updating p_z from v_z, vector length of nz * nx
*  rvdz        : sparse matrix computing dv/dz

*  RpxBvx      : FD stencil for updating p_x from v_x, vector length of nz * nx
*  rvdx        : sparse matrix computing dv/dx
"""
struct TdParams{Ti<:Int64, Tv<:AbstractFloat}

    data_format  :: DataType
    nz           :: Ti
    nx           :: Ti
    npml         :: Ti
    free_surface :: Bool
    Nz           :: Ti
    Nx           :: Ti
    ntop         :: Ti
    spt2wfd      :: Vector{Ti}
    spt2bnd      :: Vector{Ti}
    wfd2bnd      :: Vector{Ti}
    order        :: Ti
    dz           :: Tv
    dx           :: Tv
    dt           :: Tv
    tmax         :: Tv
    nt           :: Ti
    rho          :: Matrix{Tv}
    vel          :: Matrix{Tv}
    MvzBvz       :: Vector{Tv}
    MvzBp        :: Vector{Tv}
    RvzBp        :: Vector{Tv}
    dpdz         :: SparseMatrixCSC{Tv,Ti}
    rpdz         :: SparseMatrixCSC{Tv,Ti}
    MvxBvx       :: Vector{Tv}
    MvxBp        :: Vector{Tv}
    RvxBp        :: Vector{Tv}
    dpdx         :: SparseMatrixCSC{Tv,Ti}
    rpdx         :: SparseMatrixCSC{Tv,Ti}
    MpzBpz       :: Vector{Tv}
    MpzBvz       :: Vector{Tv}
    dvdz         :: SparseMatrixCSC{Tv,Ti}
    rvdz         :: SparseMatrixCSC{Tv,Ti}
    MpxBpx       :: Vector{Tv}
    MpxBvx       :: Vector{Tv}
    RpBv         :: Vector{Tv}
    dvdx         :: SparseMatrixCSC{Tv,Ti}
    rvdx         :: SparseMatrixCSC{Tv,Ti}
end

"""
   Overloading the show function for TdParams
"""
function show(io::IO, params::TdParams)

    # size of model
    @printf("nz = %4d, nx = %4d, npml = %4d\n", params.nz, params.nx, params.npml)
    @printf("dz = %4.1f, dx = %4.1f\n", params.dz, params.dx)
    @printf("time step size = %4.4f, max simulation time = %4.1f, number of step = %5d\n", params.dt, params.tmax, params.nt)

    # top boundary condition
    if params.free_surface
       @printf("top boundary conditon = %20s\n", "free surface")
    else
       @printf("top boundary conditon = %20s\n", "PML obsorbing")
    end

    # data presicion (Float32, Float64)
    @printf("data_format=%16s\n", params.data_format)

    # compute the memory cost of time domain parameters
    esize = sizeof(params)
    esize+= 172 + 8*3 + 1 + 8*3
    esize+= sizeof(params.spt2wfd) + sizeof(params.spt2bnd) + sizeof(params.wfd2bnd)
    esize+= 8*2 + sizeof(params.data_format) * 4
    esize+= sizeof(params.rho) + sizeof(params.vel)
    esize+= sizeof(params.MvzBvz)*2 + sizeof(params.MvxBvx)*2 + sizeof(params.MvzBp)*4 + sizeof(params.RvzBp)*3
    esize+= (sizeof(params.dpdz) + sizeof(params.dpdz.m)*2 + sizeof(params.dpdz.rowval) + sizeof(params.dpdz.colptr) + sizeof(params.dpdz.nzval)) * 2
    esize+= (sizeof(params.dpdx) + sizeof(params.dpdx.m)*2 + sizeof(params.dpdx.rowval) + sizeof(params.dpdx.colptr) + sizeof(params.dpdx.nzval)) * 2
    esize+= (sizeof(params.rpdz) + sizeof(params.rpdz.m)*2 + sizeof(params.rpdz.rowval) + sizeof(params.rpdz.colptr) + sizeof(params.rpdz.nzval)) * 2
    esize+= (sizeof(params.rpdz) + sizeof(params.rpdx.m)*2 + sizeof(params.rpdx.rowval) + sizeof(params.rpdx.colptr) + sizeof(params.rpdx.nzval)) * 2

    # convert to mega byte
    esize = esize / 1024^2
    @printf("esize of TdParams=%6.3f M\n", esize)

    return nothing
end

"""
   Padding PML obsorbing layers to model parameters
"""
function  model_padding(m::Matrix{Tv}, npml::Ti, free_surface::Bool) where {Ti<:Int64, Tv<:AbstractFloat}

    # model size
    (nz, nx) = size(m)

    # model size after padding
    Nx = nx + 2 * npml
    if free_surface
       Nz   = nz + npml
       ntop = 0
    else
       Nz   = nz + npml * 2
       ntop = npml
    end

    # allocate memory for the padded model
    mpad = zeros(eltype(m), Nz, Nx)

    # central part
    for ix = 1 : nx
        for iz = 1 : nz
            mpad[iz+ntop, ix+npml] = m[iz, ix]
        end
    end

    # left and right sides
    for ix = 1 : npml
        for iz = 1 : nz
            mpad[iz+ntop,ix        ] = m[iz,1 ]
            mpad[iz+ntop,ix+npml+nx] = m[iz,nx]
        end
    end

    # upper and bottom sides
    for ix = 1 : nx
        for iz = 1 : npml

            if free_surface == false              #  top
               mpad[iz,ix+npml] = m[1,ix]
            end
            mpad[iz+nz+ntop,ix+npml] = m[nz,ix]   # bottom
        end
    end

    # corners
    for ix = 1 : npml
        for iz = 1 : npml
            if free_surface == false
               mpad[iz, ix        ] = m[1, 1 ]    # upper-left
               mpad[iz, ix+npml+nx] = m[1, nx]    # upper-right
            end
            mpad[iz+ntop+nz, ix        ] = m[nz, 1 ]   # bottom-left
            mpad[iz+ntop+nz, ix+npml+nx] = m[nz, nx]   # bottom-right
        end
    end

    return mpad
end

"""
   Smoothing the 2D physical model with Hanning or Gaussian window whose half length is L
"""
function model_smooth(par::Matrix{Tv}, L::Ti; filter_flag="hanning") where {Ti<:Int64, Tv<:AbstractFloat}

    data_format = eltype(par)

    # size of smoother
    n = 2*L + 1

    # padding model before smoothing to handle boundary part
    par1 = model_padding(par, L, false)

    # convert the type of smoother
    if filter_flag == "hanning"
       s = convert(Vector{data_format}, hanning(n))
    elseif filter_flag == "gaussian"
       s = convert(Vector{data_format}, gaussian(n, 0.5))
    else
       error("wronge filter, only supper hanning and gaussian smoother")
    end
    s = s / sum(s)

    # smoothing via 2D convolution
    par1 = conv2(s,s,par1)

    # truncate to original size
    par1 = par1[2*L+1:end-2*L, 2*L+1:end-2*L]

    return par1
end

"""
compute buoyancy from density term
"""
function density2buoyancy(rho::Matrix{Tv}) where {Tv<:AbstractFloat}

    data_format = eltype(rho)

    # dimensions of the model
    (nz, nx) = size(rho)

    # allocate space for buoyancy
    buoyancy = zeros(data_format, nz, nx)

    # constant
    c_one = one(data_format)

    # compute element-wise buoyancy
    for i2 = 1 : nx
        for i1 = 1 : nz
            buoyancy[i1,i2] = c_one / rho[i1,i2]
        end
    end

    return buoyancy
end


"""
   compute lambda from velocity and density
"""
function vp2lambda(rho::Matrix{Tv}, vp::Matrix{Tv}) where {Tv<:AbstractFloat}

    (m , n ) = size(rho)
    (m1, n1) = size(vp )
    if m1 != m || n1 != n
       error("size mismatch")
    end

    lambda = zeros(eltype(vp), m, n)
    for j = 1 : n
        for i = 1 : m
            lambda[i,j] = rho[i,j] * vp[i,j] * vp[i,j]
        end
    end
    return lambda
end

"""
   get the auxillary index vectors which facilitate extracting wavefield from snapshot,
saving wavefield boundaries and insert the boundaries for wavefield reconstruction.
"""
function get_index(nz::Ti, nx::Ti, npml::Ti, free_surface::Bool; order=2) where {Ti<:Int64}

    # size of wavefield
    len1 = nz * nx

    # length of boundary values need to save
    len2 = 2*order*nz + 2*order*(nx-2*order)

    # allocate 3 vectors of integer
    index1 = zeros(Int64, len1)
    index2 = zeros(Int64, len2)
    index3 = zeros(Int64, len2)

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
            index1[icount] = tmp + i1
        end
    end

    # map snapshot to wavefield boundary
    # left side
    for i = 1 : order
        ilower = (npml+i-1) * Nz + ntop + 1
        iupper = ilower + nz - 1
        idx1   = (i-1)*nz + 1
        idx2   = idx1 +nz - 1
        index2[idx1:idx2] = collect(ilower : iupper)
    end

    # top side and bottom side
    icount = order * nz
    for ix = npml+order+1 : npml+nx-order
        tmp = (ix-1) * Nz + ntop

        for iz = 1 : order
            icount = icount + 1
            index2[icount] = tmp + iz
            icount = icount + 1
            index2[icount] = tmp + nz - iz + 1
        end
    end

    # right side
    for i = 1 : order
        ilower = (npml+nx-order+i-1) * Nz + ntop + 1
        iupper = ilower + nz - 1
        index2[icount+1: icount+nz] = collect(ilower : iupper)
        icount = icount + nz
    end

    # map wavefield to the boundary of Wavefield
    # left side
    index3[1:order*nz] = collect(1:order*nz)

    # top and bottom part in the middle
    icount = order * nz
    for ix = order+1 : nx-order
        tmp = (ix-1) * nz

        for iz = 1 : order
            icount = icount + 1
            index3[icount] = tmp + iz
            icount = icount + 1
            index3[icount] = tmp + nz - iz + 1
        end
    end

    # right side
    index3[end-order*nz+1:end] = collect((nx-order)*nz+1 : nz*nx)

    return index1, index2, index3

end

"""
   compute the finite-difference coefficients based on Tylor series expansion
"""
function finite_difference_coefficients(M::Integer)

    if M < 1
       error("M must be bigger or equal to 1")
    elseif M == 1
       c = ones(1)
       return c
    elseif M >= 2
       c = zeros(M)
       for i = 1 : M
           left  = (-1)^(i+1) / (2*i-1)
           right = 1.0
           for j = 1 : M
               if j != i
                  right = right * abs((2*j-1)^2 / ((2*i-1)^2 - (2*j-1)^2))
               end
           end
           c[i] = left * right
       end
       return c
    end
end

# finite difference coefficients estimated by minimizing LS of dispersion equation
const fdc = [[0.1129042e1, -0.4301412e-1],
             [0.1185991e1, -0.7249965e-1, 0.6301572e-2],
             [0.1217990e1, -0.9382142e-1, 0.1507536e-1, -0.1700324e-2],
             [0.1236607e1, -0.1082265   , 0.2343440e-1, -0.5033546e-2, 0.6817483e-3],
             [0.1247662e1, -0.1175538   , 0.2997970e-1, -0.8719078e-2, 0.2215897e-2, -0.3462075e-3],
             [0.1254799e1, -0.1238928   , 0.3494371e-1, -0.1208897e-1, 0.4132531e-2, -0.1197110e-2, 0.2122227e-3],
             [0.1259312e1, -0.1280347   , 0.3841945e-1, -0.1473229e-1, 0.5924913e-2, -0.2248618e-2, 0.7179226e-3, -0.1400855e-3],
             [0.1262502e1, -0.1310244   , 0.4103928e-1, -0.1686807e-1, 0.7530520e-2, -0.3345071e-2, 0.1380367e-2, -0.4808410e-3, 0.1023759e-3],
             [0.1264748e1, -0.1331606   , 0.4296909e-1, -0.1851897e-1, 0.8861071e-2, -0.4347073e-2, 0.2076101e-2, -0.9164925e-3, 0.3437446e-3, -0.7874250e-4]];

"""
   test the stable condition of finite difference method
"""
function is_stable(vel::Matrix{Tv}, dz::Tv, dx::Tv, dt::Tv) where {Tv<:AbstractFloat}

    vmax = maximum(vel)
    beta = vmax * sqrt(1/(dz^2) + 1/(dx^2))
    tt   = 1.0 / beta

    if beta * dt >= 1.0
       println("maximum acceptable time step: $tt")
       println("current dt is $dt")
       error("decrease dt or increase dz or dx")
    end

    return true
end

"""
   the constructor for TdParams
"""
function TdParams(rho, vel, free_surface::Bool, dz, dx, dt, tmax;
         data_format=Float32, fd_flag="taylor", order=2, npml=20, apml=900.) where {Ti<:Int64}

    # get the size of the model
    (nz, nx) = size(rho)
    if (nz, nx) != size(vel)
       error("model size mismatch")
    end

    # unify the data format
    rho = convert(Matrix{data_format}, rho)
    vel = convert(Matrix{data_format}, vel)
    dz  = convert(data_format, dz)
    dx  = convert(data_format, dx)
    dt  = convert(data_format, dt)
    tmax= convert(data_format, tmax)

    # test the stable condition
    is_stable(vel, dz, dx, dt)

    # compute finite difference coefficients
    if fd_flag == "taylor"
       fd_coefficients = finite_difference_coefficients(order)
    elseif fd_flag == "ls"
       fd_coefficients = fdc[order-1]
    else
       error("non-supported method for computing FD coefficients")
    end
    fd_coefficients = convert(Vector{data_format}, fd_coefficients)
    order = length(fd_coefficients)

    # size of model after padding PML layers
    if free_surface
       Nz   = nz +     npml
       ntop = 0
    else
       Nz   = nz + 2 * npml
       ntop = npml
    end
    Nx = nx + 2 * npml

    # number of time stepping
    nt = round(Int64, tmax/dt) + 1

    # compute auxillary index vectors
    (spt2wfd, spt2bnd, wfd2bnd) = get_index(nz, nx, npml, free_surface; order=order)

    # Z direction pml coefficients
    if free_surface
       (gammaz_integer, gammaz_stagger) = get_gamma_profile(Nz, npml, dz; apml=apml, oneside_flag=true)
    else
       (gammaz_integer, gammaz_stagger) = get_gamma_profile(Nz, npml, dz; apml=apml, oneside_flag=false)
    end

    # X direction pml coefficients
    (gammax_integer, gammax_stagger) = get_gamma_profile(Nx, npml, dx; apml=apml, oneside_flag=false)

    # data format convert
    gammaz_integer = convert(Vector{data_format}, gammaz_integer)
    gammaz_stagger = convert(Vector{data_format}, gammaz_stagger)
    gammax_integer = convert(Vector{data_format}, gammax_integer)
    gammax_stagger = convert(Vector{data_format}, gammax_stagger)

    # compute buoyancy and bulk modulus
    buoy0 = density2buoyancy(rho)
    bulk0 = vp2lambda(rho, vel)

    # implement free surface condition
    if free_surface
       for j = 1 : nx
           bulk0[1,j] = data_format(0.0)
       end
    end

    # padding the model parameters to accomdate PML
    buoy1 = model_padding(buoy0, npml, free_surface)
    bulk1 = model_padding(bulk0, npml, free_surface)

    # FD stencil with PML boundaries
    (MvzBvz, MvzBp , dpdz) = Mvz(gammaz_stagger, buoy1, dz, dt, fd_coefficients)
    (MvxBvx, MvxBp , dpdx) = Mvx(gammax_stagger, buoy1, dx, dt, fd_coefficients)
    (MpzBpz, MpzBvz, dvdz) = Mpz(gammaz_integer, bulk1, dz, dt, fd_coefficients)
    (MpxBpx, MpxBvx, dvdx) = Mpx(gammax_integer, bulk1, dx, dt, fd_coefficients)

    # FD stencil with rigid boundaries
    (RvzBp , rpdz) = Rvz(buoy0, dz, dt, fd_coefficients)
    (RvxBp , rpdx) = Rvx(buoy0, dx, dt, fd_coefficients)
    (RpBv, rvdz, rvdx) = Rpzx(bulk0, dz, dx, dt, fd_coefficients)

    # call the default struct constructor
    return TdParams(data_format, nz, nx, npml, free_surface, Nz, Nx, ntop,
                       spt2wfd, spt2bnd, wfd2bnd, order,
                       dz, dx, dt, tmax, nt, rho, vel,
                       MvzBvz, MvzBp , RvzBp , dpdz, rpdz,
                       MvxBvx, MvxBp , RvxBp , dpdx, rpdx,
                       MpzBpz, MpzBvz,         dvdz, rvdz,
                       MpxBpx, MpxBvx, RpBv  , dvdx, rvdx)
end


# """
#    immutable struct contain the model parameters
# *  data_format : Float32 or Float64
# *  nz          : depth of model
# *  nx          : length of model
# *  npml        : number of PML layers
# *  free_surface: true -> no obsorbtion on the surface
#                  false-> PML on the surface
# *  Nz          : depth of model (include PML padding)
# *  Nx          : length of model (include PML padding)
# *  ntop        : number of padded PML layers on the top
# *  spt2wfd     : auxillary index vector for computing wavefield from snapshot
# *  spt2bnd     : auxillary index vector for computing wavefield boundaries from snapshot
# *  bnd2wfd     : auxillary index vector for inserting wavefield boundaries to wavefield
# *  dz          : vertical grid size
# *  dx          : horizontal grid size
# *  dt          : time step size
# *  tmax        : maximum modelling length
# *  nt          : number of time steps
# *  rho         : 2D density model
# *  vel         : 2D P-wave velocity model
# *  fd_coefficients : finite difference coefficients
# *  order       : the order of precision, will decide the number of layers required for wavefield reconstruction
# """
# struct TdParams{Ti<:Int64, Tv<:AbstractFloat}
#
#     data_format  :: DataType
#     nz           :: Ti
#     nx           :: Ti
#     npml         :: Ti
#     free_surface :: Bool
#     Nz           :: Ti
#     Nx           :: Ti
#     ntop         :: Ti
#     spt2wfd      :: Vector{Ti}
#     spt2bnd      :: Vector{Ti}
#     bnd2wfd      :: Vector{Ti}
#     dz           :: Tv
#     dx           :: Tv
#     dt           :: Tv
#     tmax         :: Tv
#     nt           :: Ti
#     rho          :: Matrix{Tv}
#     vel          :: Matrix{Tv}
#     fd_coefficients :: Vector{Tv}
# end
#
# """
#    Overloading the show function for TdParams
# """
# function show(io::IO, params::TdParams)
#
#     # size of model
#     @printf("nz = %4d, nx = %4d, npml = %4d\n", params.nz, params.nx, params.npml)
#     @printf("dz = %4.1f, dx = %4.1f\n", params.dz, params.dx)
#     @printf("time step size = %4.4f, max simulation time = %4.1f, number of step = %5d\n", params.dt, params.tmax, params.nt)
#
#     # top boundary condition
#     if params.free_surface
#        @printf("top boundary conditon = %20s\n", "free surface")
#     else
#        @printf("top boundary conditon = %20s\n", "PML obsorbing")
#     end
#
#     # data presicion (Float32, Float64)
#     @printf("data_format=%16s\n", params.data_format)
#
#     return nothing
# end
#
# """
#    Padding PML obsorbing layers to model parameters
# """
# function  model_padding(m::Matrix{Tv}, npml::Ti, free_surface::Bool) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # model size
#     (nz, nx) = size(m)
#
#     # model size after padding
#     Nx = nx + 2 * npml
#     if free_surface
#        Nz   = nz + npml
#        ntop = 0
#     else
#        Nz   = nz + npml * 2
#        ntop = npml
#     end
#
#     # allocate memory for the padded model
#     mpad = zeros(eltype(m), Nz, Nx)
#
#     # central part
#     for ix = 1 : nx
#         for iz = 1 : nz
#             mpad[iz+ntop, ix+npml] = m[iz, ix]
#         end
#     end
#
#     # left and right sides
#     for ix = 1 : npml
#         for iz = 1 : nz
#             mpad[iz+ntop,ix        ] = m[iz,1 ]
#             mpad[iz+ntop,ix+npml+nx] = m[iz,nx]
#         end
#     end
#
#     # upper and bottom sides
#     for ix = 1 : nx
#         for iz = 1 : npml
#
#             if free_surface == false              #  top
#                mpad[iz,ix+npml] = m[1,ix]
#             end
#             mpad[iz+nz+ntop,ix+npml] = m[nz,ix]   # bottom
#         end
#     end
#
#     # corners
#     for ix = 1 : npml
#         for iz = 1 : npml
#             if free_surface == false
#                mpad[iz, ix        ] = m[1, 1 ]    # upper-left
#                mpad[iz, ix+npml+nx] = m[1, nx]    # upper-right
#             end
#             mpad[iz+ntop+nz, ix        ] = m[nz, 1 ]   # bottom-left
#             mpad[iz+ntop+nz, ix+npml+nx] = m[nz, nx]   # bottom-right
#         end
#     end
#
#     return mpad
# end
#
# """
#    Smoothing the 2D physical model with Hanning or Gaussian window whose half length is L
# """
# function model_smooth(par::Matrix{Tv}, L::Ti; filter_flag="hanning") where {Ti<:Int64, Tv<:AbstractFloat}
#
#     data_format = eltype(par)
#
#     # size of smoother
#     n = 2*L + 1
#
#     # padding model before smoothing to handle boundary part
#     par1 = model_padding(par, L, false)
#
#     # convert the type of smoother
#     if filter_flag == "hanning"
#        s = convert(Vector{data_format}, hanning(n))
#     elseif filter_flag == "gaussian"
#        s = convert(Vector{data_format}, gaussian(n, 0.5))
#     else
#        error("wronge filter, only supper hanning and gaussian smoother")
#     end
#     s = s / sum(s)
#
#     # smoothing via 2D convolution
#     par1 = conv2(s,s,par1)
#
#     # truncate to original size
#     par1 = par1[2*L+1:end-2*L, 2*L+1:end-2*L]
#
#     return par1
# end
#
# """
#    get the auxillary index vectors which facilitate extracting wavefield from snapshot,
# saving wavefield boundaries and insert the boundaries for wavefield reconstruction.
# """
# function get_index(nz::Ti, nx::Ti, npml::Ti, free_surface::Bool; order=2) where {Ti<:Int64}
#
#     # size of wavefield
#     len1 = nz * nx
#
#     # length of boundary values need to save
#     len2 = 2*order*nz + 2*order*(nx-2*order)
#
#     # allocate 3 vectors of integer
#     index1 = zeros(Int64, len1)
#     index2 = zeros(Int64, len2)
#     index3 = zeros(Int64, len2)
#
#     # size of padded model
#     if free_surface
#        Nz = nz + npml
#        ntop = 0
#     else
#        Nz = nz + npml * 2
#        ntop = npml
#     end
#     Nx = nx + npml * 2
#
#     # map snapshot to wavefield
#     icount = 0
#     for i2 = 1 : nx
#         tmp= (i2 + npml - 1)*Nz + ntop
#
#         for i1 = 1 : nz
#             icount = icount + 1
#             index1[icount] = tmp + i1
#         end
#     end
#
#     # map snapshot to wavefield boundary
#     # left side
#     for i = 1 : order
#         ilower = (npml+i-1) * Nz + ntop + 1
#         iupper = ilower + nz - 1
#         idx1   = (i-1)*nz + 1
#         idx2   = idx1 +nz - 1
#         index2[idx1:idx2] = collect(ilower : iupper)
#     end
#
#     # top side and bottom side
#     icount = order * nz
#     for ix = npml+order+1 : npml+nx-order
#         tmp = (ix-1) * Nz + ntop
#
#         for iz = 1 : order
#             icount = icount + 1
#             index2[icount] = tmp + iz
#             icount = icount + 1
#             index2[icount] = tmp + nz - iz + 1
#         end
#     end
#
#     # right side
#     for i = 1 : order
#         ilower = (npml+nx-order+i-1) * Nz + ntop + 1
#         iupper = ilower + nz - 1
#         index2[icount+1: icount+nz] = collect(ilower : iupper)
#         icount = icount + nz
#     end
#
#     # map wavefield to the boundary of Wavefield
#     # left side
#     index3[1:order*nz] = collect(1:order*nz)
#
#     # top and bottom part in the middle
#     icount = order * nz
#     for ix = order+1 : nx-order
#         tmp = (ix-1) * nz
#
#         for iz = 1 : order
#             icount = icount + 1
#             index3[icount] = tmp + iz
#             icount = icount + 1
#             index3[icount] = tmp + nz - iz + 1
#         end
#     end
#
#     # right side
#     index3[end-order*nz+1:end] = collect((nx-order)*nz+1 : nz*nx)
#
#     return index1, index2, index3
#
# end
#
# """
#    compute the finite-difference coefficients based on Tylor series expansion
# """
# function finite_difference_coefficients(M::Integer)
#
#     if M < 1
#        error("M must be bigger or equal to 1")
#     elseif M == 1
#        c = ones(1)
#        return c
#     elseif M >= 2
#        c = zeros(M)
#        for i = 1 : M
#            left  = (-1)^(i+1) / (2*i-1)
#            right = 1.0
#            for j = 1 : M
#                if j != i
#                   right = right * abs((2*j-1)^2 / ((2*i-1)^2 - (2*j-1)^2))
#                end
#            end
#            c[i] = left * right
#        end
#        return c
#     end
# end
#
# # finite difference coefficients estimated by minimizing LS of dispersion equation
# const fdc = [[0.1129042e1, -0.4301412e-1],
#              [0.1185991e1, -0.7249965e-1, 0.6301572e-2],
#              [0.1217990e1, -0.9382142e-1, 0.1507536e-1, -0.1700324e-2],
#              [0.1236607e1, -0.1082265   , 0.2343440e-1, -0.5033546e-2, 0.6817483e-3],
#              [0.1247662e1, -0.1175538   , 0.2997970e-1, -0.8719078e-2, 0.2215897e-2, -0.3462075e-3],
#              [0.1254799e1, -0.1238928   , 0.3494371e-1, -0.1208897e-1, 0.4132531e-2, -0.1197110e-2, 0.2122227e-3],
#              [0.1259312e1, -0.1280347   , 0.3841945e-1, -0.1473229e-1, 0.5924913e-2, -0.2248618e-2, 0.7179226e-3, -0.1400855e-3],
#              [0.1262502e1, -0.1310244   , 0.4103928e-1, -0.1686807e-1, 0.7530520e-2, -0.3345071e-2, 0.1380367e-2, -0.4808410e-3, 0.1023759e-3],
#              [0.1264748e1, -0.1331606   , 0.4296909e-1, -0.1851897e-1, 0.8861071e-2, -0.4347073e-2, 0.2076101e-2, -0.9164925e-3, 0.3437446e-3, -0.7874250e-4]];
#
# """
#    test the stable condition of finite difference method
# """
# function is_stable(vel::Matrix{Tv}, dz::Tv, dx::Tv, dt::Tv) where {Tv<:AbstractFloat}
#
#     vmax = maximum(vel)
#     beta = vmax * sqrt(1/(dz^2) + 1/(dx^2))
#     tt   = 1.0 / beta
#
#     if beta * dt >= 1.0
#        println("maximum acceptable time step: $tt")
#        println("current dt is $dt")
#        error("decrease dt or increase dz or dx")
#     end
#
#     return true
# end
#
# """
#    the constructor for TdParams
# """
# function TdParams(rho, vel, npml::Ti, free_surface::Bool, dz, dx, dt, tmax;
#          data_format=Float32, fd_flag="taylor", order=2) where {Ti<:Int64}
#
#     # get the size of the model
#     (nz, nx) = size(rho)
#     if (nz, nx) != size(vel)
#        error("model size mismatch")
#     end
#
#     # unify the data format
#     rho = convert(Matrix{data_format}, rho)
#     vel = convert(Matrix{data_format}, vel)
#     dz  = convert(data_format, dz)
#     dx  = convert(data_format, dx)
#     dt  = convert(data_format, dt)
#     tmax= convert(data_format, tmax)
#
#     # test the stable condition
#     is_stable(vel, dz, dx, dt)
#
#     # compute finite difference coefficients
#     if fd_flag == "taylor"
#        fd_coefficients = finite_difference_coefficients(order)
#     elseif fd_flag == "ls"
#        fd_coefficients = fdc[order-1]
#     else
#        error("non-supported method for computing FD coefficients")
#     end
#     fd_coefficients = convert(Vector{data_format}, fd_coefficients)
#     order = length(fd_coefficients)
#
#     # size of model after padding PML layers
#     if free_surface
#        Nz   = nz +     npml
#        ntop = 0
#     else
#        Nz   = nz + 2 * npml
#        ntop = npml
#     end
#     Nx = nx + 2 * npml
#
#     # number of time stepping
#     nt = round(Int64, tmax/dt) + 1
#
#     # compute auxillary index vectors
#     (spt2wfd, spt2bnd, bnd2wfd) = get_index(nz, nx, npml, free_surface; order=order)
#
#     # call the default struct constructor
#     return TdParams(data_format, nz, nx, npml, free_surface, Nz, Nx, ntop,
#            spt2wfd, spt2bnd, bnd2wfd, dz, dx, dt, tmax, nt, rho, vel, fd_coefficients)
#
# end
