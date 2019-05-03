"""
   Define the partial derivative operator for v_z with PML boundary
"""
function Mvz(gammaz_stagger::Vector{Tv}, buoy::Matrix{Tv},
         dz::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(buoy)
    if nz != length(gammaz_stagger)
       error("size mismatch")
    end

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(buoy)
    double      = data_format(2.0)

    # the coefficients in front of v_z
    MvzBvz = ones(data_format, nz)
    for i = 1 : nz
        MvzBvz[i] = (double - dt * gammaz_stagger[i]) / (double + dt * gammaz_stagger[i])
    end

    # the coefficients in front of p
    MvzBp  = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1

            # averge of density
            if i < nz
               bavg = (buoy[i,j] + buoy[i+1,j]) / double
            else
               bavg = buoy[i,j]
            end

            MvzBp[idx] = (double * bavg * dt) / (double+dt*gammaz_stagger[i])
        end
    end

    # spatial partial derivative
    dpdz = zeros(data_format, nz, nz)
    for ir = 1 : nz
        lower = ir - order + 1 >= 1  ? ir-order+1 : 1
        upper = ir + order     <= nz ? ir+order   : nz

        for ic = lower : upper
            if ic-ir <= 0
               idx = ir-ic + 1
               dpdz[ir,ic] = - fd[idx] / dz
            elseif ic-ir >= 1
               idx = ic-ir
               dpdz[ir,ic] =   fd[idx] / dz
            end
        end
    end

    return MvzBvz, MvzBp, sparse(dpdz)
end

"""
   Define the partial derivative operator for v_z with rigid boundary
"""
function Rvz(buoy::Matrix{Tv}, dz::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(buoy)

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(buoy)
    double      = data_format(2.0)

    # the coefficients in front of p
    RvzBp  = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1

            # averge of density
            if i < nz
               bavg = (buoy[i,j] + buoy[i+1,j]) / double
            else
               bavg = buoy[i,j]
            end

            RvzBp[idx] = bavg * dt
        end
    end

    # spatial partial derivative
    rpdz = zeros(data_format, nz, nz)
    for ir = 1 : nz
        lower = ir - order + 1 >= 1  ? ir-order+1 : 1
        upper = ir + order     <= nz ? ir+order   : nz

        for ic = lower : upper
            if ic-ir <= 0
               idx = ir-ic + 1
               rpdz[ir,ic] = - fd[idx] / dz
            elseif ic-ir >= 1
               idx = ic-ir
               rpdz[ir,ic] =   fd[idx] / dz
            end
        end
    end

    return RvzBp, sparse(rpdz)
end

"""
   Define the partial derivative operator for v_x with PML boundary
"""
function Mvx(gammax_stagger::Vector{Tv}, buoy::Matrix{Tv},
         dx::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(buoy)
    if nx != length(gammax_stagger)
       error("size mismatch")
    end

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(buoy)
    double      = data_format(2.0)

    # the coefficients in front of v_x
    MvxBvx = ones(data_format, nx)
    for j = 1 : nx
        MvxBvx[j] = (double - dt * gammax_stagger[j]) / (double + dt * gammax_stagger[j])
    end

    # the coefficients in front of p
    MvxBp  = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1

            # averge of density
            if j < nx
               bavg = (buoy[i,j] + buoy[i,j+1]) / double
            else
               bavg = buoy[i,j]
            end

            MvxBp[idx] = (double * bavg * dt) / (double + dt * gammax_stagger[j])
        end
    end

    # spatial partial derivative dpdx
    dpdx = zeros(data_format, nx, nx)
    for ir = 1 : nx
        lower = ir - order + 1 >= 1  ? ir-order+1 : 1
        upper = ir + order     <= nx ? ir+order   : nx

        for ic = lower : upper
            if ic-ir <= 0
               idx = ir - ic + 1
               dpdx[ir,ic] = - fd[idx] / dx
            elseif ic-ir >= 1
               idx = ic - ir
               dpdx[ir,ic] =   fd[idx] / dx
            end
        end
    end

    return MvxBvx, MvxBp, sparse(dpdx)
end

"""
   Define the partial derivative operator for v_x with rigid boundary
"""
function Rvx(buoy::Matrix{Tv}, dx::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(buoy)

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(buoy)
    double      = data_format(2.0)

    # the coefficients in front of p
    RvxBp  = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1

            # averge of density
            if j < nx
               bavg = (buoy[i,j] + buoy[i,j+1]) / double
            else
               bavg = buoy[i,j]
            end

            RvxBp[idx] = bavg * dt
        end
    end

    # spatial partial derivative dpdx
    rpdx = zeros(data_format, nx, nx)
    for ir = 1 : nx
        lower = ir - order + 1 >= 1  ? ir-order+1 : 1
        upper = ir + order     <= nx ? ir+order   : nx

        for ic = lower : upper
            if ic-ir <= 0
               idx = ir - ic + 1
               rpdx[ir,ic] = - fd[idx] / dx
            elseif ic-ir >= 1
               idx = ic - ir
               rpdx[ir,ic] =   fd[idx] / dx
            end
        end
    end

    return RvxBp, sparse(rpdx)
end

"""
   Define the partial derivative operator for p_z with PML boundary
"""
function Mpz(gammaz_integer::Vector{Tv}, bulk::Matrix{Tv},
         dz::Tv, dt::Tv, fd::Vector) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(bulk)
    if nz != length(gammaz_integer)
       error("size mismatch")
    end

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(bulk)
    double      = data_format(2.0)

    # the coefficients in front of p_z
    MpzBpz = ones(data_format, nz)
    for i = 1 : nz
        MpzBpz[i] = (double - dt * gammaz_integer[i]) / (double + dt * gammaz_integer[i])
    end

    # the coefficients in front of v_z
    MpzBvz = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1
            MpzBvz[idx] = (double * bulk[i,j] * dt) / (double + dt * gammaz_integer[i])
        end
    end

    # spatial partial derivative dvdz
    dvdz = zeros(data_format, nz, nz)
    for ir = 1 : nz
        lower = ir - order     >= 1  ? ir-order   : 1
        upper = ir + order - 1 <= nz ? ir+order-1 : nz

        for ic = lower : upper
            if ic-ir <= -1
               idx = ir-ic
               dvdz[ir,ic] = - fd[idx] / dz
            elseif ic-ir >= 0
               idx = ic-ir+1
               dvdz[ir,ic] =   fd[idx] / dz
            end
        end
    end

    return MpzBpz, MpzBvz, sparse(dvdz)
end

"""
   Define the partial derivative operator for p_x with PML boundary
"""
function Mpx(gammax_integer::Vector{Tv}, bulk::Matrix{Tv},
         dx::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(bulk)
    if nx != length(gammax_integer)
       error("size mismatch")
    end

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(bulk)
    double      = data_format(2.0)

    # the coefficients in front of p_x
    MpxBpx = ones(data_format, nx)
    for j = 1 : nx
        MpxBpx[j] = (double - dt * gammax_integer[j]) / (double + dt * gammax_integer[j])
    end

    # the coefficients in front of v_x
    MpxBvx = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1
            MpxBvx[idx] = (double * bulk[i,j] * dt) / (double + dt * gammax_integer[j])
        end
    end

    # spatial partial derivative dvdx
    dvdx = zeros(data_format, nx, nx)
    for ir = 1 : nx
        lower = ir - order     >= 1  ? ir-order   : 1
        upper = ir + order - 1 <= nx ? ir+order-1 : nx

        for ic = lower : upper
            if ic-ir <= -1
               idx = ir - ic
               dvdx[ir,ic] = - fd[idx] / dx
            elseif ic-ir >= 0
               idx = ic - ir + 1
               dvdx[ir,ic] =   fd[idx] / dx
            end
        end
    end

    return MpxBpx, MpxBvx, sparse(dvdx)
end

"""
   Define the partial derivative operator for p_z with rigid boundary
"""
function Rpzx(bulk::Matrix{Tv}, dz::Tv, dx::Tv, dt::Tv, fd::Vector) where {Tv<:AbstractFloat}

    # size of the model
    (nz, nx) = size(bulk)

    # precision order
    order = length(fd)

    # data format
    data_format = eltype(bulk)
    double      = data_format(2.0)

    # the coefficients in front of v_z
    RpBv = ones(data_format, nz*nx)
    idx = 0
    for j = 1 : nx
        for i = 1 : nz
            idx = idx + 1
            RpBv[idx] = bulk[i,j] * dt
        end
    end

    # spatial partial derivative dvdz
    rvdz = zeros(data_format, nz, nz)
    for ir = 1 : nz
        lower = ir - order     >= 1  ? ir-order   : 1
        upper = ir + order - 1 <= nz ? ir+order-1 : nz

        for ic = lower : upper
            if ic-ir <= -1
               idx = ir-ic
               rvdz[ir,ic] = - fd[idx] / dz
            elseif ic-ir >= 0
               idx = ic-ir+1
               rvdz[ir,ic] =   fd[idx] / dz
            end
        end
    end

    # spatial partial derivative dvdx
    rvdx = zeros(data_format, nx, nx)
    for ir = 1 : nx
        lower = ir - order     >= 1  ? ir-order   : 1
        upper = ir + order - 1 <= nx ? ir+order-1 : nx

        for ic = lower : upper
            if ic-ir <= -1
               idx = ir - ic
               rvdx[ir,ic] = - fd[idx] / dx
            elseif ic-ir >= 0
               idx = ic - ir + 1
               rvdx[ir,ic] =   fd[idx] / dx
            end
        end
    end

    return RpBv, sparse(rvdz), sparse(rvdx)
end

# """
#    sparse matrix for finite-difference grid with rigid boundary
# """
# struct RigidFDStencil{Ti<:Int64, Tv<:AbstractFloat}
#     MvzBp  :: SparseMatrixCSC{Tv,Ti}
#     MvxBp  :: SparseMatrixCSC{Tv,Ti}
#     MpzBvz :: SparseMatrixCSC{Tv,Ti}
#     MpxBvx :: SparseMatrixCSC{Tv,Ti}
# end
#
# """
#    constructor for finite difference stencil with rigid boundary condition
# """
# function RigidFDStencil(params::ModelParams)
#
#     # expand physical model parameters
#     lambda = vp2lambda(params.vel, params.rho)
#
#     MvzBp  = Mvz(params.rho, params.dz, params.dt, params.fd_coefficients)
#     MvxBp  = Mvx(params.rho, params.dx, params.dt, params.fd_coefficients)
#     MpzBvz = Mpz(lambda    , params.dz, params.dt, params.fd_coefficients)
#     MpxBvx = Mpx(lambda    , params.dx, params.dt, params.fd_coefficients)
#
#     return RigidFDStencil(MvzBp, MvxBp, MpzBvz, MpxBvx)
# end


# """
#    sparse matrix for finite-difference grid with PML boundary
# """
# struct ObsorbFDStencil{Ti<:Int64, Tv<:AbstractFloat}
#     MvzBvz :: Vector{Tv}
#     MvzBp  :: SparseMatrixCSC{Tv,Ti}
#     MvxBvx :: Vector{Tv}
#     MvxBp  :: SparseMatrixCSC{Tv,Ti}
#     MpzBpz :: Vector{Tv}
#     MpzBvz :: SparseMatrixCSC{Tv,Ti}
#     MpxBpx :: Vector{Tv}
#     MpxBvx :: SparseMatrixCSC{Tv,Ti}
# end
#
# """
#    sparse matrix for finite-difference grid with rigid boundary
# """
# struct RigidFDStencil{Ti<:Int64, Tv<:AbstractFloat}
#     MvzBp  :: SparseMatrixCSC{Tv,Ti}
#     MvxBp  :: SparseMatrixCSC{Tv,Ti}
#     MpzBvz :: SparseMatrixCSC{Tv,Ti}
#     MpxBvx :: SparseMatrixCSC{Tv,Ti}
# end
#
# """
#    compute lambda from velocity and density
# """
# function vp2lambda(vp::Matrix{Tv}, rho::Matrix{Tv}) where {Tv<:AbstractFloat}
#
#     (m , n ) = size(vp)
#     (m1, n1) = size(rho)
#     if m1 != m || n1 != n
#        error("size mismatch")
#     end
#
#     lambda = zeros(eltype(vp), m, n)
#     for j = 1 : n
#         for i = 1 : m
#             lambda[i,j] = rho[i,j] * vp[i,j]^2
#         end
#     end
#     return lambda
# end
#
# """
#    Define the partial derivative operator for v_z with PML boundary
# """
# function Mvz(vz_pml::SparseMatrixCSC{Tv,Ti}, rho::Matrix{Tv},
#          dz::Tv, dt::Tv, fd::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # data format
#     data_format = eltype(rho)
#     double_dt   = data_format(2.0 * dt)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dz
#
#     # model size
#     (m, n) = size(vz_pml)
#     c_diag = zeros(data_format, m*n)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#
#             if iz < m
#                a = (rho[iz+1,ix]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vz_pml[iz,ix] / 2
#             denum[(ix-1)*m+iz] = 1.0   / (a+b)
#             c_diag[(ix-1)*m+iz]= (a-b) / (a+b)
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, m, m)
#     for iz = 1 : m
#         lower = iz - order + 1 >= 1 ? iz-order+1 : 1
#         upper = iz + order     <= m ? iz+order   : m
#
#         for ix = lower : upper
#             if ix-iz <= 0
#                idx = iz-ix + 1
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 1
#                idx = ix-iz
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MvzBvz = c_diag
#     MvzBp  = spdiagm(0 => denum) * kron(spdiagm(0=>ones(data_format,n)), sparse(tmp))
#
#     return MvzBvz, MvzBp
#
# end
#
# """
#    Define the partial derivative operator for v_z with rigid boundary
# """
# function Mvz(rho::Matrix{Tv}, dz::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     # data format
#     data_format = eltype(rho)
#     double_dt   = data_format(2.0 * dt)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dz
#
#     # model size
#     (m, n) = size(rho)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#
#             if iz < m
#                a = (rho[iz+1,ix]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             denum[(ix-1)*m+iz] = 1.0 / a
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, m, m)
#     for iz = 1 : m
#         lower = iz - order + 1 >= 1 ? iz-order+1 : 1
#         upper = iz + order     <= m ? iz+order   : m
#
#         for ix = lower : upper
#             if ix-iz <= 0
#                idx = iz-ix + 1
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 1
#                idx = ix-iz
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MvzBp  = spdiagm(0 => denum) * kron(spdiagm(0=>ones(data_format,n)), sparse(tmp))
#     return MvzBp
# end
#
# """
#    Define the partial derivative operator for v_x with PML boundary
# """
# function Mvx(vx_pml::SparseMatrixCSC{Tv,Ti}, rho::Matrix{Tv},
#          dx::Tv, dt::Tv, fd::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # data format
#     data_format  = eltype(rho)
#     double_dt    = data_format(2.0 * dt)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dx
#
#     # model size
#     (m, n) = size(vx_pml)
#     c_diag = zeros(data_format, m*n)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#
#             if ix < n
#                a = (rho[iz,ix+1]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vx_pml[iz,ix] / 2
#             denum[(ix-1)*m+iz] = 1     / (a+b)
#             c_diag[(ix-1)*m+iz]= (a-b) / (a+b)
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, n, n)
#     for iz = 1 : n
#         lower = iz - order + 1 >= 1 ? iz-order+1 : 1
#         upper = iz + order     <= n ? iz+order   : n
#
#         for ix = lower : upper
#             if ix-iz <= 0
#                idx = iz-ix + 1
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 1
#                idx = ix-iz
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MvxBvx = c_diag
#     MvxBp  = spdiagm(0 => denum) * kron(sparse(tmp), spdiagm(0=>ones(data_format,m)))
#
#     return MvxBvx, MvxBp
#
# end
#
# """
#    Define the partial derivative operator for v_x with rigid boundary
# """
# function Mvx(rho::Matrix{Tv}, dx::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     # data format
#     data_format  = eltype(rho)
#     double_dt = data_format(2.0 * dt)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dx
#
#     # model size
#     (m, n) = size(rho)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             # average density
#             if ix < n
#                a = (rho[iz,ix+1]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             denum[(ix-1)*m+iz] = 1.0 / a
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, n, n)
#     for iz = 1 : n
#         lower = iz - order + 1 >= 1 ? iz-order+1 : 1
#         upper = iz + order     <= n ? iz+order   : n
#
#         for ix = lower : upper
#             if ix-iz <= 0
#                idx = iz-ix + 1
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 1
#                idx = ix-iz
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MvxBp  = spdiagm(0 => denum) * kron(sparse(tmp), spdiagm(0=>ones(data_format,m)))
#     return MvxBp
# end
#
# """
#    Define the partial derivative operator for p_z with PML boundary
# """
# function Mpz(pz_pml::SparseMatrixCSC{Tv,Ti}, lambda::Matrix{Tv},
#          dz::Tv, dt::Tv, fd::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # data format
#     data_format = eltype(lambda)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dz
#
#     # model size
#     (m, n) = size(pz_pml)
#     c_diag = zeros(data_format, m*n)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + pz_pml[iz,ix]/2
#             b = 1/dt - pz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             c_diag[(ix-1)*m+iz]= b / a
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, m, m)
#     for iz = 1 : m
#         lower = iz - order     >= 1 ? iz-order   : 1
#         upper = iz + order - 1 <= m ? iz+order-1 : m
#
#         for ix = lower : upper
#             if ix-iz <= -1
#                idx = iz-ix
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 0
#                idx = ix-iz+1
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MpzBpz = c_diag
#     MpzBvz = spdiagm(0 => denum) * kron(spdiagm(0=>ones(data_format,n)), sparse(tmp))
#
#     return MpzBpz, MpzBvz
# end
#
# """
#    Define the partial derivative operator for p_z with rigid boundary
# """
# function Mpz(lambda::Matrix{Tv}, dz::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     # data format
#     data_format = eltype(lambda)
#     c           = fd / dz
#
#     # model size
#     (m, n) = size(lambda)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             denum[(ix-1)*m+iz] = lambda[iz,ix] * dt
#         end
#     end
#
#     # spatial partial derivative
#     order = length(fd)
#     tmp   = zeros(data_format, m, m)
#     for iz = 1 : m
#         lower = iz - order     >= 1 ? iz-order   : 1
#         upper = iz + order - 1 <= m ? iz+order-1 : m
#
#         for ix = lower : upper
#             if ix-iz <= -1
#                idx = iz-ix
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 0
#                idx = ix-iz+1
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MpzBvz = spdiagm(0 => denum) * kron(spdiagm(0=>ones(data_format,n)), sparse(tmp))
#     return MpzBvz
#
# end
#
# """
#    Define the partial derivative operator for p_x with PML boundary
# """
# function Mpx(px_pml::SparseMatrixCSC{Tv,Ti}, lambda::Matrix{Tv},
#          dx::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat, Ti<:Int64}
#
#     # data format
#     data_format = eltype(lambda)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dx
#
#     (m, n) = size(px_pml)
#     c_diag = zeros(data_format, m*n)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + px_pml[iz,ix]/2
#             b = 1/dt - px_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             c_diag[(ix-1)*m+iz]= b / a
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, n, n)
#     for iz = 1 : n
#         lower = iz - order     >= 1 ? iz-order   : 1
#         upper = iz + order - 1 <= n ? iz+order-1 : n
#
#         for ix = lower : upper
#             if ix-iz <= -1
#                idx = iz-ix
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 0
#                idx = ix-iz+1
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MpxBpx = c_diag
#     MpxBvx = spdiagm(0 => denum) * kron(sparse(tmp), spdiagm(0=>ones(data_format,m)))
#
#     return  MpxBpx, MpxBvx
# end
#
# """
#    Define the partial derivative operator for p_x with rigid boundary
# """
# function Mpx(lambda::Matrix{Tv}, dx::Tv, dt::Tv, fd::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     # data format
#     data_format = eltype(lambda)
#
#     # finite-difference coefficients
#     order = length(fd)
#     c     = fd / dx
#
#     (m, n) = size(lambda)
#     denum  = zeros(data_format, m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             denum[(ix-1)*m+iz] = lambda[iz,ix] * dt
#         end
#     end
#
#     # spatial partial derivative
#     tmp = zeros(data_format, n, n)
#     for iz = 1 : n
#         lower = iz - order     >= 1 ? iz-order   : 1
#         upper = iz + order - 1 <= n ? iz+order-1 : n
#
#         for ix = lower : upper
#             if ix-iz <= -1
#                idx = iz-ix
#                tmp[iz,ix] = - c[idx]
#             elseif ix-iz >= 0
#                idx = ix-iz+1
#                tmp[iz,ix] =   c[idx]
#             end
#         end
#     end
#
#     MpxBvx = spdiagm(0 => denum) * kron(sparse(tmp), spdiagm(0=>ones(data_format,m)))
#
#     return MpxBvx
# end
#
# """
#    constructor for obsorbing finite difference stencil
# """
# function ObsorbFDStencil(params::ModelParams)
#
#     # pml coefficient
#     pml = PMLCoefficients(params)
#
#     # padding physical model parameters
#     lambda = vp2lambda(params.vel, params.rho)
#     lambda = model_padding(lambda    , params.npml, params.free_surface)
#     rho    = model_padding(params.rho, params.npml, params.free_surface)
#
#     (MvzBvz, MvzBp ) = Mvz(pml.vz, rho   , params.dz, params.dt, params.fd_coefficients)
#     (MvxBvx, MvxBp ) = Mvx(pml.vx, rho   , params.dx, params.dt, params.fd_coefficients)
#     (MpzBpz, MpzBvz) = Mpz(pml.pz, lambda, params.dz, params.dt, params.fd_coefficients)
#     (MpxBpx, MpxBvx) = Mpx(pml.px, lambda, params.dx, params.dt, params.fd_coefficients)
#
#     return ObsorbFDStencil(MvzBvz, MvzBp , MvxBvx, MvxBp,
#                            MpzBpz, MpzBvz, MpxBpx, MpxBvx)
# end
#
# """
#    constructor for finite difference stencil with rigid boundary condition
# """
# function RigidFDStencil(params::ModelParams)
#
#     # expand physical model parameters
#     lambda = vp2lambda(params.vel, params.rho)
#
#     MvzBp  = Mvz(params.rho, params.dz, params.dt, params.fd_coefficients)
#     MvxBp  = Mvx(params.rho, params.dx, params.dt, params.fd_coefficients)
#     MpzBvz = Mpz(lambda    , params.dz, params.dt, params.fd_coefficients)
#     MpxBvx = Mpx(lambda    , params.dx, params.dt, params.fd_coefficients)
#
#     return RigidFDStencil(MvzBp, MvxBp, MpzBvz, MpxBvx)
# end

# old version
# function compute_relative_error(fd::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     M = length(fd)
#
#     interval = 0.001 * pi
#     beta = collect(interval : interval : pi)
#     N    = length(beta)
#     err  = zeros(N)
#
#     for i = 1 : N
#         tmp = 0.0
#         for m = 1 : M
#             tmp = tmp + fd[m] * sin((m-0.5)*beta[i])
#         end
#         err[i] = 2.0 * tmp / beta[i] - 1.0
#     end
#
#     return beta, err
# end


# the coefficients obtained by minimize absolute error
# fdc = Vector{Vector{Float64}}(9)
# fdc[1] = [0.1129136e1, -0.4304542e-1]
# fdc[2] = [0.1186247e1, -0.7266808e-1, 0.6351497e-2]
# fdc[3] = [0.1218159e1, -0.9397218e-1, 0.1519043e-1, -0.1742128e-2]
# fdc[4] = [0.1236425e1, -0.1081130   , 0.2339911e-1, -0.5061550e-2, 0.7054313e-3]
# fdc[5] = [0.1247576e1, -0.1174969   , 0.2997288e-1, -0.8741572e-2, 0.2262285e-2, -0.3745306e-3]
# fdc[6] = [0.1254380e1, -0.1235307   , 0.3467231e-1, -0.1192915e-1, 0.4057090e-2, -0.1191005e-2, 0.2263204e-3]
# fdc[7] = [0.1259012e1, -0.1277647   , 0.3820715e-1, -0.1458251e-1, 0.5845385e-2, -0.2213861e-2, 0.7243880e-3, -0.1566173e-3]
# fdc[8] = [0.1262147e1, -0.1306967   , 0.4075792e-1, -0.1665221e-1, 0.7377057e-2, -0.3258150e-2, 0.1336259e-2, -0.4775830e-3, 0.1151664e-3]
# fdc[9] = [0.1264362e1, -0.1327958   , 0.4264687e-1, -0.1824918e-1, 0.8656223e-2, -0.4200034e-2, 0.1989180e-2, -0.8686637e-3, 0.3342741e-3, -0.8854090e-4]

# M = 5
#
# (x_axis, err ) = compute_relative_error(fdc[M-1])
# (x_axis, err1) = compute_relative_error(rfdc[M-1])
#
# cfdc = finite_difference_coefficient(M)
# (x_axis, err2) = compute_relative_error(cfdc)
#
# plot(x_axis, err, label="new_a")
# plot(x_axis, err1, label="new_r")
# plot(x_axis, err2, label="old"); legend();

# """
#    finite difference stencil for updating vz
# """
# function Mvz(vz_pml::SparseMatrixCSC{Tv,Ti}, rho::Matrix{Tv},
#          dz::Tv, dt::Tv) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     data_format  = eltype(rho)
#     (m, n) = size(vz_pml)
#     double_dt = data_format(2.0 * dt)
#
#     a1 = data_format( 9 / 8 )
#     a2 = data_format(-1 / 24)
#     c1 = a1 / dz
#     c2 = a2 / dz
#
#     c3    = zeros(data_format, m*n)
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             if iz < m
#                a = (rho[iz+1,ix]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vz_pml[iz,ix] / 2
#
#             denum[(ix-1)*m+iz] = 1.0   / (a+b)
#             c3[(ix-1)*m+iz]    = (a-b) / (a+b)
#         end
#     end
#
#     tmp = zeros(data_format, m, m)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#
#     for iz = 2: m-2
#         tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
#         tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
#     end
#
#     tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
#     tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
#
#     MvzBvz = c3
#     MvzBp  = spdiagm(denum) * kron(speye(data_format, n), sparse(tmp))
#
#     return MvzBvz, MvzBp
#
# end
#
# """
#    rigid boundary
# """
# function Mvz(rho::Matrix{Tv}, dz::Tv, dt::Tv) where {Tv<:AbstractFloat}
#
#     data_format  = eltype(rho)
#     (m, n) = size(rho)
#     double_dt = data_format(2.0 * dt)
#
#     a1 = data_format( 9 / 8 )
#     a2 = data_format(-1 / 24)
#     c1 = a1 / dz
#     c2 = a2 / dz
#
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             if iz < m
#                a = (rho[iz+1,ix]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             denum[(ix-1)*m+iz] = 1 / a
#         end
#     end
#
#     tmp = zeros(data_format, m, m)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#
#     for iz = 2: m-2
#         tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
#         tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
#     end
#
#     tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
#     tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
#
#     MvzBp  = spdiagm(denum) * kron(speye(data_format, n), sparse(tmp))
#
#     return MvzBp
#
# end

# """
#    finite difference stencil for updating vx
# """
# function Mvx(vx_pml::SparseMatrixCSC{Tv,Ti}, rho::Matrix{Tv},
#          dx::Tv, dt::Tv) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     data_format  = eltype(rho)
#     (m, n) = size(vx_pml)
#     double_dt = data_format(2.0 * dt)
#
#     a1 = data_format(9/8  )
#     a2 = data_format(-1/24)
#     c1 = a1 / dx
#     c2 = a2 / dx
#
#     c3    = zeros(data_format, m*n)
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             if ix < n
#                a = (rho[iz,ix+1]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vx_pml[iz,ix] / 2
#
#             denum[(ix-1)*m+iz] = 1     / (a+b)
#             c3[(ix-1)*m+iz]    = (a-b) / (a+b)
#         end
#     end
#
#     tmp = zeros(data_format, n, n)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#
#     for ix = 2: n-2
#         tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
#         tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
#     end
#
#     tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
#     tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
#
#     MvxBvx = c3
#     MvxBp = spdiagm(denum) * kron(sparse(tmp), speye(data_format, m))
#
#     return MvxBvx, MvxBp
#
# end
#
# """
#    rigid boundary Mvx
# """
# function Mvx(rho::Matrix{Tv}, dx::Tv, dt::Tv) where {Tv<:AbstractFloat}
#
#     data_format  = eltype(rho)
#     (m, n) = size(rho)
#     double_dt = data_format(2.0 * dt)
#
#     a1 = data_format(9/8  )
#     a2 = data_format(-1/24)
#     c1 = a1 / dx
#     c2 = a2 / dx
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             if ix < n
#                a = (rho[iz,ix+1]+rho[iz,ix]) / double_dt
#             else
#                a = rho[iz,ix] / dt
#             end
#             denum[(ix-1)*m+iz] = 1 / a
#         end
#     end
#
#     tmp = zeros(data_format, n, n)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#
#     for ix = 2: n-2
#         tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
#         tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
#     end
#
#     tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
#     tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
#
#     MvxBp = spdiagm(denum) * kron(sparse(tmp), speye(data_format, m))
#
#     return MvxBp
#
# end

# """
#    finite difference stencil for updating pz
# """
# function Mpz(pz_pml::SparseMatrixCSC{Tv,Ti}, lambda::Matrix{Tv},
#          dz::Tv, dt::Tv) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     data_format  = eltype(lambda)
#     (m, n) = size(pz_pml)
#
#     a1 = data_format( 9/8 )
#     a2 = data_format(-1/24)
#     c1 = a1 / dz
#     c2 = a2 / dz
#
#     c3    = zeros(data_format, m*n)
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + pz_pml[iz,ix]/2
#             b = 1/dt - pz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             c3[(ix-1)*m+iz]= b / a
#         end
#     end
#
#     tmp = zeros(data_format, m, m)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#
#     for iz = 3: m-1
#         tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
#         tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
#     end
#
#     tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
#
#     MpzBpz = c3
#     MpzBvz = spdiagm(denum) * kron(speye(data_format, n), sparse(tmp))
#
#     return MpzBpz, MpzBvz
#
# end

# """
#    rigid boundary Mpz with high-order finite difference approximation
# """
# function Mpz(lambda::Matrix{Tv}, dz::Tv, dt::Tv) where {Tv<:AbstractFloat}
#
#     data_format  = eltype(lambda)
#     (m, n) = size(lambda)
#
#     a1 = data_format( 9/8 )
#     a2 = data_format(-1/24)
#     c1 = a1 / dz
#     c2 = a2 / dz
#
#     c3    = zeros(data_format, m*n)
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             denum[(ix-1)*m+iz] = lambda[iz,ix] * dt
#         end
#     end
#
#     tmp = zeros(data_format, m, m)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#
#     for iz = 3: m-1
#         tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
#         tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
#     end
#
#     tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
#
#     MpzBvz = spdiagm(denum) * kron(speye(data_format, n), sparse(tmp))
#
#     return MpzBvz
# end

# """
#    finite difference stencil for updating px
# """
# function Mpx(px_pml::SparseMatrixCSC{Tv,Ti}, lambda::Matrix{Tv},
#          dx::Tv, dt::Tv) where {Tv<:AbstractFloat, Ti<:Int64}
#
#     data_format  = eltype(lambda)
#     (m, n) = size(px_pml)
#
#     a1 = data_format( 9/8 )
#     a2 = data_format(-1/24)
#     c1 = a1 / dx
#     c2 = a2 / dx
#
#     c3    = zeros(data_format, m*n)
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + px_pml[iz,ix]/2
#             b = 1/dt - px_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             c3[(ix-1)*m+iz]= b / a
#         end
#     end
#
#     tmp = zeros(data_format, n, n)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#
#     for ix = 3: n-1
#         tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
#         tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
#     end
#
#     tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
#
#     MpxBpx = c3
#     MpxBvx = spdiagm(denum) * kron(sparse(tmp), speye(data_format, m))
#
#     return  MpxBpx, MpxBvx
# end
#
#
# """
#    rigid boundary Mpx
# """
# function Mpx(lambda::Matrix{Tv}, dx::Tv, dt::Tv) where {Tv<:AbstractFloat}
#
#     data_format  = eltype(lambda)
#     (m, n) = size(lambda)
#
#     a1 = data_format( 9/8 )
#     a2 = data_format(-1/24)
#     c1 = a1 / dx
#     c2 = a2 / dx
#
#     denum = zeros(data_format, m*n)
#
#     for ix = 1 : n
#         for iz = 1 : m
#             denum[(ix-1)*m+iz] = lambda[iz,ix] * dt
#         end
#     end
#
#     tmp = zeros(data_format, n, n)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#
#     for ix = 3: n-1
#         tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
#         tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
#     end
#
#     tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
#
#     MpxBvx = spdiagm(denum) * kron(sparse(tmp), speye(data_format, m))
#
#     return MpxBvx
# end

# """
#    constructor for finite difference sparse matrix
# """
# function FiniteDiffMatrix(par::PhysicalModel)
#
#     # pml coefficient
#     pml = PMLCoefficient(par)
#
#     # expand physical model parameters
#     lambda = vp2lambda(par.vel, par.rho)
#     lambda = model_expand(lambda , par.npml, par.free_surface)
#     rho    = model_expand(par.rho, par.npml, par.free_surface)
#
#     (MvzBvz, MvzBp ) = Mvz(pml.vz, rho   , par.dz, par.dt)
#     (MvxBvx, MvxBp ) = Mvx(pml.vx, rho   , par.dx, par.dt)
#     (MpzBpz, MpzBvz) = Mpz(pml.pz, lambda, par.dz, par.dt)
#     (MpxBpx, MpxBvx) = Mpx(pml.px, lambda, par.dx, par.dt)
#
#     return FiniteDiffMatrix(MvzBvz, MvzBp , MvxBvx, MvxBp,
#                             MpzBpz, MpzBvz, MpxBpx, MpxBvx)
#
# end
#
# """
#    finite difference sparse matrix with rigit boundary condition
# """
# function RigidFiniteDiffMatrix(par::PhysicalModel)
#
#     # expand physical model parameters
#     lambda = vp2lambda(par.vel, par.rho)
#
#     MvzBp  = Mvz(par.rho, par.dz, par.dt)
#     MvxBp  = Mvx(par.rho, par.dx, par.dt)
#     MpzBvz = Mpz(lambda , par.dz, par.dt)
#     MpxBvx = Mpx(lambda , par.dx, par.dt)
#
#     return RigidFiniteDiffMatrix(MvzBp , MvxBp,
#                                  MpzBvz, MpxBvx)
#
# end
