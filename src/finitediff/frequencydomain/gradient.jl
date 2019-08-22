"""
   compute the gradient of objective function with respect to velocity model for one
frequency slice, r is the residue obtained by subtracting observed data from synthetic,
u is source side wavefield with PML boundary part, H is the LU factorization
of Helmhotz operator.
"""
function velocity_gradient(r::MonochromaticRecordings, u::Vector{Tv},
         H, params::FdParams) where {Tv<:Complex{Float64}}

    # insert the conplex conjugate of residue to snapshot
    b = Vector{Complex{Float64}}(undef, params.Nz*params.Nx)
    inject_rec2spt!(b, r)

    # square of radian frequency
    omega2 = r.omega * r.omega

    # some coefficients of mixed-grid
    w1 = 0.6287326
    w2 = 0.0928166
    w3 = (1.0 - w1 - 4.0*w2) / 4.0

    # create a table save the coefficients for mass acceleration average
    C  = OffsetArray{Float64}(undef, -1:1, -1:1)

    C[-1,-1] = w3; C[-1,0] = w2; C[-1,1] = w3;
    C[ 0,-1] = w2; C[ 0,0] = w1; C[ 0,1] = w2;
    C[ 1,-1] = w3; C[ 1,0] = w2; C[ 1,1] = w3;

    # compute adjoint wavefield
    u_adj = zeros(Complex{Float64}, params.Nz*params.Nx)
    ldiv!(u_adj, transpose(H), b)

    # apply the imaging condition to get the gradient
    g = zeros(Complex{Float64}, params.nz * params.nx)

    # both u and u_adj are vectors include PML layers
    for i2 = 1 : params.nx
        idx2 = (i2-1) * params.nz

        for i1 = 1 : params.nz
            idx1 = idx2 + i1

            # the constant part
            tmp   = params.spt2wfd[idx1]
            alpha = -2.0 * omega2 * u[tmp]
            alpha = alpha / (params.rho[i1,i2] * (params.vel[i1,i2])^3)

            for d2 = -1 : 1
                l2 = (i2+params.npml+d2-1) * params.Nz

                for d1 = -1 : 1
                    if 1 <= i1+d1 || !params.free_surface
                       l1 = l2 + i1 + d1 + params.ntop

                       g[idx1] = g[idx1] + C[d1,d2] * u_adj[l1]
                    end
                end
            end

            g[idx1] = g[idx1] * alpha
        end
    end

    return real(g)
end

# some other test the transpose of LU decomposition
# A = sprand(Complex{Float64}, 10, 10, 0.3) + spdiagm(0=>rand(Complex{Float64}, 10));
# B = Matrix(transpose(A))
# B = Sparse(B
# H = lu(B)
# P = lu(A)
#
# x0 = rand(Complex{Float64}, 10);
# b  = Vector{Complex{Float64}}(undef, 10);
# mul!(b, transpose(A), x0);
#
# x = Vector{Complex{Float64}}(undef, 10);
# ldiv!(x, H, b)
#
# x1 = Vector{Complex{Float64}}(undef, 10);
# ldiv!(x1, transpose(P), b)
