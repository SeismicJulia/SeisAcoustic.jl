function CGLS(A::Function, b::Ts; x0="NULL", path="NULL", params=[],
         shift=0.0, tol=1e-6, maxIter=50, print_flag=false) where {Ts<:Union{String, Vector{String}}}

    if path == "NULL"
       path = pwd()
    end

    # create folder to save files of model size
    path_model = joinpath(path, "model")
    if !isdir(path_model)
       mkdir(path_model)
       if !isdir(path_model)
          error("could not create directory")
       end
    end
    x = join([path_model "/x.bin"])
    s = join([path_model "/s.bin"])
    p = join([path_model "/p.bin"])

    # create folders to save residue
    path_residue   = joinpath(path, "residue")
    if !isdir(path_residue)
       mkdir(path_residue)
       if !isdir(path_residue)
          error("could not create directory")
       end
    end

    # create folders to save synthetic data with current model
    path_forward = joinpath(path, "forward")
    if !isdir(path_forward)
       mkdir(path_forward)
       if !isdir(path_forward)
          error("could not create directory")
       end
    end

    # intermidiate variables
    if typeof(b) == String
       r = join([path_residue "/res.bin"])
       q = join([path_forward "/fwd.bin"])
    elseif typeof(b) == Vector{String}
       n = length(b)
       r = Vector{String}(n)
       q = Vector{String}(n)
       for i = 1 : n
           r[i] = join([path_residue "/res" "_" "$i" ".bin"])
           q[i] = join([path_forward "/fwd" "_" "$i" ".bin"])
       end
    end

    # provided initial guess
    if x0 != "NULL"

       cp(x0, x, remove_destination=true)
       dcal = A(x, 1; params=params)
       x_plus_alpha_q!(r, b, -1.0, dcal)

       if typeof(params) <: Dict{Symbol, Any}
          params[:path_m] = s
       elseif typeof(params) == Vector{Dict}
         for i = 1 : length(params)
             params[i][:path_m] = s
         end
       end
       s = A(r, 2; params=params)
       x_plus_alpha_q!(s, -shift, x)

    # zero vector as initial guess
    else

       cp(b, r, remove_destination=true)
       if typeof(params) <: Dict{Symbol, Any}
          params[:path_m] = s
       elseif typeof(params) == Vector{Dict}
          for i = 1 : length(params)
              params[i][:path_m] = s
          end
       end
       s = A(r, 2; params=params)

       # create x
       (hdr, ds) = read_USdata(s)
       tmp = zeros(ds)
       write_USdata(x, hdr, tmp)
    end

    # compute residue
    data_fitting = (norm(r))^2
    constraint   = 0.0
    cost0 = data_fitting + constraint
    convergence = Float64[]; push!(convergence, 1.0);

    # initialize some intermidiate vectors
    cp(s, p, remove_destination=true)

    norms0= norm(s)
    gamma = norms0^2
    normx = norm(x)
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

    if typeof(params) <: Dict{Symbol, Any}
       params[:path_fwd] = q
    elseif typeof(params) == Vector{Dict}
       for i = 1 : length(params)
           params[i][:path_fwd] = q[i]
       end
    end

    while (k < maxIter) && run_flag

          k = k + 1

          q = A(p, 1; params=params)

          delta = (norm(q))^2 + shift * (norm(p))^2
          indefinite = delta <= 0 ? true : false
          delta      = delta == 0.? eps(): delta

          alpha = gamma / delta

          x_plus_alpha_q!(x,  alpha, p)
          x_plus_alpha_q!(r, -alpha, q)

          data_fitting = (norm(r))^2
          constraint   = shift * (norm(x))^2
          cost = data_fitting + constraint

          # save the intermidiate result
          path_iter = join([path_model "/iteration" "_" "$k" ".bin"])
          cp(x, path_iter, remove_destination=true)

          s = A(r, 2; params=params)
          x_plus_alpha_q!(s, -shift, x)

          norms  = norm(s)
          gamma0 = gamma
          gamma  = norms^2
          beta   = gamma / gamma0

          (hdr, ds) = read_USdata(s)
          (hdr, dp) = read_USdata(p)
          dp = ds + beta * dp
          write_USdata(p, hdr, dp)

          # check the stopping crietia
          normx = norm(x)
          xmax  = normx > xmax ? normx : xmax
          if norms <= norms0 * tol || normx * tol >= 1.0
             run_flag = false
          end

          # print information
          resNE = cost / cost0
          if print_flag
             @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k, data_fitting, constraint, normx, resNE);
          end
          push!(convergence, resNE);
    end

    return x, convergence

end


function build_model_window(iz::Ti, hl::Ti, phi::PhysicalModel) where (Ti<:Int64)

    if hl == 0
       return ones(phi.data_format, phi.nz, phi.nx)
    else
       if iz < hl + 1
          error("too small iz")
       end

       window = zeros(phi.data_format, phi.nz)
       taper  = hanning(2*hl+1)
       window[iz-hl:iz] = taper[1:hl+1]

       scale = repmat(window, 1, phi.nx)
       scale[iz+1:end,:] = one(phi.data_format)

       return scale
    end
end

"""
   wrapper for forward born approximationa and write the recordings in standard
uniform-sampled data format.
"""
function wrap_preconditioned_born_forward(par::Dict)

    # current model parameter
    (hdr, delta_lambda) = read_USdata(par[:path_m])
    delta_lambda = vec(delta_lambda)

    # apply the preconditioner
    (hdr, scale) = read_USdata(par[:path_pre])
    for i = 1 : length(delta_lambda)
                         #pesudo hessian model windowing
        delta_lambda[i] = scale[i] * par[:model_window][i] * delta_lambda[i]
    end

    # boundary of source-side wavefield and the last one
    (bnd, wfd) = read_bound_wavefield(par[:path_bnd])

    # Born forward modeling
    rec = born_approximation_forward(par[:irz], par[:irx], delta_lambda, par[:fidMtx], par[:src], bnd, par[:fidMtxT], par[:phi])

    # write the recording to hard drive
    write_record(par[:path_fwd], rec, par[:phi])
end

"""
  wrapper for the adjoint of born approximationa and write the image in standard
uniform-sampled data format.
"""
function wrap_preconditioned_born_adjoint(par::Dict)

    # receiver location
    rec = read_record(par[:path_obs], par[:phi], 350, 50)

    # boundary of source-side wavefield and the last one
    (bnd, wfd) = read_bound_wavefield(par[:path_bnd])

    # adjoint of Born approximation
    delta_lambda = born_approximation_adjoint(rec, par[:fidMtx], par[:src], bnd, wfd, par[:fidMtxT], par[:phi])

    # apply the adjoint of preconditioner
    (hdr, scale) = read_USdata(par[:path_pre])
    for i = 1 : length(delta_lambda)
                         #model windowing       pesudo hessian
        delta_lambda[i] = par[:model_window][i] * scale[i] * delta_lambda[i]
    end
    # write the recording to hard drive
    write_image(par[:path_adj], delta_lambda, par[:phi])

end

function wrap_preconditioned_born_approximation(x::Ts, iflag::Ti; params=[]) where {Ts<:Union{String, Vector{String}}, Ti<:Int64}

    # number of shot
    ns = length(params)

    # forward modeling operator
    if iflag == 1

       # change the path of current model
       for i = 1 : ns
           params[i][:path_m] = x
       end

       # the result are write to par[i][:path_fwd]
       pmap(wrap_preconditioned_born_forward, params)

       # return the path
       d = Vector{String}(ns)
       for i = 1 : ns
           d[i] = params[i][:path_fwd]
       end
       return d

    # adjoint operator
    elseif iflag == 2

       for i = 1 : ns
           params[i][:path_obs] = x[i]
       end
       pmap(wrap_preconditioned_born_adjoint, params)
       stack_image(params)

       return params[1][:path_m]
    end

end
