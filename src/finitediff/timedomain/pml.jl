"""
   struct for the coefficients of the split PML boundary
"""
struct PMLCoefficients{Tv<:AbstractFloat, Ti<:Int64}
    vz::SparseMatrixCSC{Tv, Ti}
    vx::SparseMatrixCSC{Tv, Ti}
    pz::SparseMatrixCSC{Tv, Ti}
    px::SparseMatrixCSC{Tv, Ti}
end

"""
   constructor for PML obsorbing boundary coefficients
"""
function PMLCoefficients(params::ModelParams)

    data_format = params.data_format

    # coefficient depends on number of layers
    R = data_format(0.01)
    if params.npml == 10
       R = data_format(0.001)
    elseif params.npml >= 20
       R = data_format(0.0001)
    end

    ll = data_format(2.)
    vmin = minimum(params.vel)
    alpha_z = data_format(-1.5) * vmin / (params.npml * params.dz) * log(R)
    alpha_x = data_format(-1.5) * vmin / (params.npml * params.dx) * log(R)

    # allocate memory for PML coefficient
    vz = zeros(data_format, params.Nz, params.Nx)
    vx = zeros(data_format, params.Nz, params.Nx)
    pz = zeros(data_format, params.Nz, params.Nx)
    px = zeros(data_format, params.Nz, params.Nx)

    # expand model to accomdate PML layers
    rho = model_padding(params.rho, params.npml, params.free_surface)

    if params.free_surface
       #  lower
       for iz = params.nz+1 : params.Nz
           vz[iz,:] = alpha_z * ((iz+0.5-params.nz) / params.npml)^ll * rho[iz,:]
           pz[iz,:] = alpha_z * ((iz    -params.nz) / params.npml)^ll * ones(data_format, params.Nx)
       end

       #  left
       for ix = 1 : params.npml
           vx[:,ix] = alpha_x * ((params.npml+0.5-ix) / params.npml)^ll * rho[:,ix]
           px[:,ix] = alpha_x * ((params.npml+1  -ix) / params.npml)^ll * ones(data_format, params.Nz)
       end

       #  right
       for ix = params.nx + params.npml + 1 : params.Nx
           vx[:,ix] = alpha_x * ((ix+0.5-params.nx-params.npml) / params.npml)^ll * rho[:,ix]
           px[:,ix] = alpha_x * ((ix    -params.nx-params.npml) / params.npml)^ll * ones(data_format, params.Nz)
       end
    else
       #  upper and left
       for i = 1 : params.npml
           vz[i,:] = alpha_z * ((params.npml+0.5-i) / params.npml)^ll * rho[i,:]
           vx[:,i] = alpha_x * ((params.npml+0.5-i) / params.npml)^ll * rho[:,i]
           pz[i,:] = alpha_z * ((params.npml+1  -i) / params.npml)^ll * ones(data_format, params.Nx)
           px[:,i] = alpha_x * ((params.npml+1  -i) / params.npml)^ll * ones(data_format, params.Nz)
       end

       #  lower
       for iz = params.nz + params.npml + 1 : params.Nz
           vz[iz,:] = alpha_z * ((iz+0.5-params.nz-params.npml) / params.npml)^ll * rho[iz,:]
           pz[iz,:] = alpha_z * ((iz    -params.nz-params.npml) / params.npml)^ll * ones(data_format, params.Nx)
       end

       #  right
       for ix = params.nx + params.npml + 1 : params.Nx
           vx[:,ix] = alpha_x * ((ix+0.5-params.nx-params.npml) / params.npml)^ll * rho[:,ix]
           px[:,ix] = alpha_x * ((ix    -params.nx-params.npml) / params.npml)^ll * ones(data_format, params.Nz)
       end
    end

    return PMLCoefficients(sparse(vz), sparse(vx), sparse(pz), sparse(px))
end
