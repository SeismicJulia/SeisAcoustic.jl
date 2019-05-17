"""
   set up damping profile for PML boundary conditions, n is the model size include pml,
npml is the number of PML layers, omega is radian frequency.
"""
function get_gamma_profile(n::Ti, npml::Ti, h; apml=900.0, oneside_flag=false) where {Ti<:Int64}

    # constants
    pi_over_2 = 0.5 * pi   # pi/2

    # damping profile of PML
    gamma_integer = zeros(Float64, n)
    gamma_stagger = zeros(Float64, n)

    # PML width
    pml_width = npml * h

    # model width include PML layers
    model_width = (n - 1) * h

    # both side of gamma_integer and one side of gamma_stagger
    for i = 1 : npml

        # current coordinate of x
        x_integer = (i-1) * h
        x_stagger = (i-1+0.5) * h

        # distance of current grid to the edge of computation area
        dis_integer = pml_width - x_integer
        dis_stagger = pml_width - x_stagger

        # pml coefficients
        gamma_integer[i] = apml * (1.0 - cos(dis_integer / pml_width * pi_over_2))
        gamma_stagger[i] = apml * (1.0 - cos(dis_stagger / pml_width * pi_over_2))

        # symmetrical gamma_integer
        gamma_integer[n-i+1] = gamma_integer[i]
    end

    # one side PML obsorbing boundary
    if oneside_flag
       for i = 1 : npml
           gamma_integer[i] = 0.0
           gamma_stagger[i] = 0.0
       end
    end

    # the right side of gamma_stagger
    for i = 1 : npml + 1

        # the coordinate of the stagger grid (from right to left)
        x_stagger = model_width + 0.5 * h - (i-1)*h

        # distance of current grid to the edge of computation area
        dis_stagger = x_stagger - (model_width - pml_width)

        # pml coefficients
        gamma_stagger[n-i+1] = apml * (1.0 - cos(dis_stagger / pml_width * pi_over_2))
    end

    return gamma_integer, gamma_stagger
end

# """
#    struct for the coefficients of the split PML boundary
# """
# struct PMLCoefficients{Tv<:AbstractFloat, Ti<:Int64}
#     vz::SparseMatrixCSC{Tv, Ti}
#     vx::SparseMatrixCSC{Tv, Ti}
#     pz::SparseMatrixCSC{Tv, Ti}
#     px::SparseMatrixCSC{Tv, Ti}
# end
#
# """
#    constructor for PML obsorbing boundary coefficients
# """
# function PMLCoefficients(params::TdParams)
#
#     data_format = params.data_format
#
#     # coefficient depends on number of layers
#     R = data_format(0.01)
#     if params.npml == 10
#        R = data_format(0.001)
#     elseif params.npml >= 20
#        R = data_format(0.0001)
#     end
#
#     ll = data_format(2.)
#     vmin = minimum(params.vel)
#     alpha_z = data_format(-1.5) * vmin / (params.npml * params.dz) * log(R)
#     alpha_x = data_format(-1.5) * vmin / (params.npml * params.dx) * log(R)
#
#     # allocate memory for PML coefficient
#     vz = zeros(data_format, params.Nz, params.Nx)
#     vx = zeros(data_format, params.Nz, params.Nx)
#     pz = zeros(data_format, params.Nz, params.Nx)
#     px = zeros(data_format, params.Nz, params.Nx)
#
#     # expand model to accomdate PML layers
#     rho = model_padding(params.rho, params.npml, params.free_surface)
#
#     if params.free_surface
#        #  lower
#        for iz = params.nz+1 : params.Nz
#            vz[iz,:] = alpha_z * ((iz+0.5-params.nz) / params.npml)^ll * rho[iz,:]
#            pz[iz,:] = alpha_z * ((iz    -params.nz) / params.npml)^ll * ones(data_format, params.Nx)
#        end
#
#        #  left
#        for ix = 1 : params.npml
#            vx[:,ix] = alpha_x * ((params.npml+0.5-ix) / params.npml)^ll * rho[:,ix]
#            px[:,ix] = alpha_x * ((params.npml+1  -ix) / params.npml)^ll * ones(data_format, params.Nz)
#        end
#
#        #  right
#        for ix = params.nx + params.npml + 1 : params.Nx
#            vx[:,ix] = alpha_x * ((ix+0.5-params.nx-params.npml) / params.npml)^ll * rho[:,ix]
#            px[:,ix] = alpha_x * ((ix    -params.nx-params.npml) / params.npml)^ll * ones(data_format, params.Nz)
#        end
#     else
#        #  upper and left
#        for i = 1 : params.npml
#            vz[i,:] = alpha_z * ((params.npml+0.5-i) / params.npml)^ll * rho[i,:]
#            vx[:,i] = alpha_x * ((params.npml+0.5-i) / params.npml)^ll * rho[:,i]
#            pz[i,:] = alpha_z * ((params.npml+1  -i) / params.npml)^ll * ones(data_format, params.Nx)
#            px[:,i] = alpha_x * ((params.npml+1  -i) / params.npml)^ll * ones(data_format, params.Nz)
#        end
#
#        #  lower
#        for iz = params.nz + params.npml + 1 : params.Nz
#            vz[iz,:] = alpha_z * ((iz+0.5-params.nz-params.npml) / params.npml)^ll * rho[iz,:]
#            pz[iz,:] = alpha_z * ((iz    -params.nz-params.npml) / params.npml)^ll * ones(data_format, params.Nx)
#        end
#
#        #  right
#        for ix = params.nx + params.npml + 1 : params.Nx
#            vx[:,ix] = alpha_x * ((ix+0.5-params.nx-params.npml) / params.npml)^ll * rho[:,ix]
#            px[:,ix] = alpha_x * ((ix    -params.nx-params.npml) / params.npml)^ll * ones(data_format, params.Nz)
#        end
#     end
#
#     return PMLCoefficients(sparse(vz), sparse(vx), sparse(pz), sparse(px))
# end
