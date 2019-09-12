function matrix_mul(y, x, iflag::Ti; A::Matrix)

    # forward
    if iflag == 1
       mul!(y, A, x)

    # adjoint
    elseif iflag == 2
       mul!(y, transpose(A), x)

    else
       error("undefined operation")
    end

    return y
end

# set up the problem
for T in (Float32, Float64, ComplexF32, ComplexF64)
end

m = 100; n = 100;
A  = rand(T, m, n)
x0 = rand(T, n)
d  = A * x0

(x, his) = cgls(matrix_mul!, y; op_params=params, maxIter=m)


"""
   CGLS algorithm for linear operators whose input and output can be save in memory
"""
function cgls(A::Tf, b::Vector{Tv}; op_params::NamedTuple=NamedTuple(),
              x0=[], shift=0.0, tol=1e-6, maxIter=50, print_flag=true) where {Tf<:Union{Function, Matrix}, Tv<:Number}

    # determine the type of A(matrix or operator)
    if isa(A, Matrix)
       explicitA = true

    elseif isa(A, Function)
       explicitA = false

    else
       error("A can be either a matrix or a function")
    end

    # initial guess of the unknowns
    if explicitA

       (m, n) = size(A)

       # provided initial guess
       if x0 != 0.0
          x = copy(x0)
          r = b - A * x
          s = A' * r - shift * x
       # zero vector as initial guess
       else
          x = zeros(eltype(b), n)
          r = copy(b)
          s = A' * r
       end

    else

       m = length(b)

       # provided initial guess
       if x0 != 0.0
          x = copy(x0)
          r = b - A(x, 1; op_params...)
          s = A(r, 2; op_params...) - shift * x
          n = length(s)
       # zero vector as initial guess
       else
          r = copy(b)
          s = A(b, 2; params...)
          n = length(s)
          x = zeros(eltype(b), n)
       end

    end

    data_fitting = dot(r,r)
    constraint   = 0.0
    cost0 = data_fitting + constraint

    # initialize some intermidiate vectors
    p = copy(s)
    norms0= vecnorm(s)
    gamma = norms0^2
    normx = vecnorm(x)
    xmax  = normx     # keep the initial one
    resNE = 1.0

    gamma0= copy(gamma)
    delta = 0.0

    # iteration counter and stop condition
    k = 0
    run_flag = true

    if print_flag
       header = "  k         data_fitting           constraint             normx                resNE"
       println(""); println(header);
       @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k, data_fitting, constraint, normx, resNE);
    end

    while (k < maxIter) && run_flag

          k = k + 1

          if explicitA
             q = A * p
          else
             q = A(p, 1; params...)
          end

          delta = (vecnorm(q))^2 + shift * (vecnorm(p))^2
          indefinite = delta <= 0 ? true : false
          delta      = delta == 0.? eps(): delta

          alpha = gamma / delta

          x = x + alpha * p
          r = r - alpha * q

          data_fitting = dot(r,r)
          constraint   = shift * dot(x, x)
          cost = data_fitting + constraint

          if explicitA
             s = A' * r - shift * x
          else
             s = A(r, 2; params...) - shift * x
          end

          norms  = norm(s)
          gamma0 = copy(gamma)
          gamma  = norms^2
          beta   = gamma / gamma0

          p = s + beta * p

          # check the stopping crietia
          normx = norm(x)
          xmax  = normx > xmax ? normx : xmax
          if norms <= norms0 * tol || normx * tol >= 1.0
             run_flag = false
          end

          # print information
          # resNE = norms / norms0
          resNE = cost  / cost0
          if print_flag
             @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k, data_fitting, constraint, normx, resNE);
          end

    end

    return x, resNE

end
