"""
   compute the gradient of objective function with respect to velocity model for one
frequency slice, u is source side wavefield, dres is the residue, H is the LU factorization
of Helmhotz operator
"""
function velocity_gradient(b::Vector{Tv}, u::Vector{Tv},
         H, params::FdParams) where {Tv<:Complex{Float64}}

    # create a table save the coefficients for mass acceleration average
    T = OffsetArray{Float64}(undef, -1:1, -1:1)
    w1= 0.6287326 * omega2
    w2= 0.0928166 * omega2
    w3 = (1.0 - w1 - 4.0*w2) / 4.0 * omega2
    T[-1,-1]=w3; T[-1, 0]=w2; T[-1, 1]=w3;
    T[ 0,-1]=w2; T[ 0, 0]=w1; T[ 0, 1]=w2;
    T[ 1,-1]=w3; T[ 1, 0]=w2; T[ 1, 1]=w3;

    # compute adjoint wavefield
    u_adj = zeros(Complex{Float64}, params.Nz, params.Nx)
    ldiv!(u_adj, transpose(H), b)

    # apply the imaging condition to get the gradient
    g = zeros(Complex{Float64}, params.nz * params.nx)

    # both u and u_adj are vectors include PML layers
    for i2 = 1 : params.nx
        tmp= (i2-1) * params.nz

        for i1 = 1 : params.nz
            idx= tmp + i1
            j  = params.spt2wfd[idx]
            cst= -2.0 * u[j] / params.rho[i1,i2] / (params.vel[i1,i2])^3

            for d2 = -1 : 1
                if 1 <= i2+d2 <= params.nx
                   tmp1 = (i2+d2-1) * params.nz

                   for d1 = -1 : 1
                       if 1 <= i1+d1 <= params.nz
                          idx1 = tmp1 + i1 + d1
                          j1   = params.spt2wfd[idx1]
                          g[idx] = g[idx] + T[d1,d2] * u_adj[j1]
                       end
                   end
                end
            end

            g[idx] = g[idx] * cst
        end
    end

    return real(g)

end


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

# LU decomposition of Helmhotz operator
omega = 2.0 * pi * 10.0;
H = get_helmholtz_LU(fparams, omega)

# generate observations of one frequency slice







# inject the residues
b = zeros(Complex{Float64}, params.Nz * params.Nx)
for i = 1 : length(dres.spt2rec)
    j = dres.spt2rec[i]
    b[j] = conj(dres.p[i])
end
