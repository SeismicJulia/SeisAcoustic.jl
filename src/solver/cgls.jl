function cgls(op::Function, b::Ts; op_params::NamedTuple, data_axpby::Function, model_axpby::Function,
              data_norm::Function, model_norm::Function, x0="NULL", dir_work="NULL",
              shift=0.0, tol=1e-6, maxIter=50, print_flag=false) where {Ts<:Union{String, Vector{String}}}

    # determine the working directory
    if dir_work == "NULL"
       dir_work =  pwd()
    end

    # create directory for variables in model space
    dir_model = joinpath(dir_work, "model_space")
    rm(dir_model, force=true, recursive=true)
    mkdir(dir_model)
    if !isdir(dir_model)
       error("could not create directory for variables in model space")
    end
    x = joinpath(dir_model, "x.rsf")
    s = joinpath(dir_model, "s.rsf")
    p = joinpath(dir_model, "p.rsf")

    # create directory for variables in data space
    dir_data = joinpath(dir_work, "data_space")
    rm(dir_data, force=true, recursive=true)
    mkdir(dir_data)
    if !isdir(dir_data)
       error("could not create directory for variables in data space")
    end

    # create directory for residue
    dir_residue = joinpath(dir_data, "residue")
    mkdir(dir_residue)
    if !isdir(dir_residue)
       error("could not create directory of residue")
    end

    # create directory for the synthetic with current model
    dir_forward = joinpath(dir_data, "forward")
    mkdir(dir_forward)
    if !isdir(dir_forward)
       error("could not create directory of forward synthetic")
    end

    # intermidiate variables
    r = dir_residue
    q = dir_forward

    # initial guess provided
    if x0 != "NULL"

       cp(x0, x, force=true)
       op(r, x, 1; op_params...)
       data_axpby(1.0, b, -1.0, r)      # r = b - r

       op(s, r, 2; op_params...)
       model_axpby(-shift, x, 1.0, s)   # s = s - shift*x

    # no initial guess
    else

       cp(b, r, force=true)
       op(s, r, 2; op_params...)

       # create x
       hdr = read_RSheader(s)
       write_RSdata(x, hdr, zeros(hdr))
    end

    # compute residue
    data_fitting = (data_norm(r))^2
    constraint   = 0.0
    cost0        = data_fitting + constraint
    convergence  = Float64[]; push!(convergence, 1.0);

    # initialize some intermidiate vectors
    cp(s, p, force=true)

    norms0= model_norm(s)
    gamma = norms0^2
    normx = model_norm(x)
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
          op(q, p, 1; op_params...)

          delta = (data_norm(q))^2 + shift * (model_norm(p))^2
          indefinite = delta <= 0.0 ? true  : false
          delta      = delta == 0.0 ? eps() : delta

          alpha = gamma / delta

          model_axpby(alpha, p, 1.0, x)  # x = x + alpha * p
          data_axpby(-alpha, q, 1.0, r)  # r = r - alpha * q

          data_fitting = (data_norm(r))^2
          constraint   = shift * (model_norm(x))^2
          cost         = data_fitting + constraint

          # save the intermidiate result
          path_iter = join([dir_model "/iteration" "_" "$k" ".rsf"])
          cp(x, path_iter, force=true)

          op(s, r, 2; op_params...)
          model_axpby(-shift, x, 1.0, s) # s = s - shift * x

          norms  = model_norm(s)
          gamma0 = gamma
          gamma  = norms^2
          beta   = gamma / gamma0

          model_axpby(1.0, s, beta, p)   # p = s + beta * p

          # check the stopping crietia
          normx = model_norm(x)
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
