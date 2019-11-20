# ==============================================================================
function time_shift(d::Matrix{Tv}, dither::Vector, dt) where {Tv<:Real}

    # to eliminate the influence of periodic property of DFT
    # maybe worth to padding zeros to both ends, no implemented at here.

    # get dimension of data
    (n1, n2) = size(d)
    n2 == length(dither) || error("check the length of dither")

    # fourier transform along time axis
    D = fft(d, 1)

    # radian frequency interval
    delta_w = 2 * pi * (1.0 / dt / n1)

    # compute half of the radian frequency axis
    # if n1 is even, nyquist frequency is computed, otherwise not
    nw = floor(Int64, n1/2) + 1
    w  = zeros(Tv, nw)
    for i = 1 : nw
        w[i] = (i-1) * delta_w
    end

    # apply time shift to each trace
    for i2 = 1 : n2
        for i1 = 2 : nw  # start from first non-zero frequency component

            # the conjugate property of real series
            j1 = n1 - i1 + 2
            D[i1,i2] = D[i1,i2] * exp(-im * w[i1] * dither[i2])
            D[j1,i2] = conj(D[i1,i2])
        end
    end

    # transform back to time domain
    return real(ifft(D,1))

end

function blend_forward(d1, d2, time_delay, dt)
    # apply time delay to far boat
    d3 = time_shift(d2, time_delay, dt)
    return d1 + d3
end

function blend_adjoint(d_blend, time_delay, dt)

    d1 = copy(d_blend)
    d2 = time_shift(d_blend, -time_delay, dt)

    return d1, d2
end

function hard_threshold!(d::Matrix{Complex{Tv}}, gamma::Tv) where {Tv<:AbstractFloat}

    (n1, n2) = size(d)
    for i2 = 1 : n2
        for i1 = 1 : n1
            if abs(d[i1,i2]) < gamma
               d[i1,i2] = 0.0im
            end
        end
    end
    return nothing
end

function deblending_by_inversion(d_blend::Matrix{Tv}, time_delay::Vector{Tv}, dt::Tv;
         max_iter=200, min_value=1e-5) where {Tv<:AbstractFloat}

     (n1, n2) =  size(d_blend)
     n2       == length(time_delay) || error("errors in length of time delay")

     # cost of objective function
     cost0 = norm(d_blend)

     # hard-thresholding value decay exponentially
     b = min_value^(1/max_iter)

     # initialize model parameter
     m1 = zeros(Tv, n1, n2)
     m2 = zeros(Tv, n1, n2)

     # iteration
     for iter = 1 : max_iter

         # compute residue
         r = blend_forward(m1, m2, time_delay, dt) - d_blend

         # verbose
         tmp = norm(r) / cost0
         println("iteration: $iter, relative cost: $tmp")

         # compute gradient
         (g1, g2) = blend_adjoint(r, time_delay, dt)

         # update model parameter
         m1 .= m1 .- g1
         m2 .= m2 .- g2
         m2[1:1350,:] .= 0.0

         # projection
         fraction = b^iter
         F1 = fft(m1)
         F2 = fft(m2)
         gamma1 = maximum(abs.(F1)) * fraction
         gamma2 = maximum(abs.(F2)) * fraction
         hard_threshold!(F1, gamma1)
         hard_threshold!(F2, gamma2)

         # transform back to time domain
         m1 = real(ifft(F1))
         m2 = real(ifft(F2))
     end

     return m1, m2
end

# read the data
dir_work = joinpath(homedir(), "Desktop/Data_Beijing_Workshop/Marmousi2");
path_near= joinpath(dir_work, "near_offset1_4ms.rsf");
path_far = joinpath(dir_work, "far_offset1_4ms.rsf");

(hdr, d1) = read_RSdata(path_near); d1 = convert(Matrix{Float64}, d1);
(hdr, d2) = read_RSdata(path_far);  d2 = convert(Matrix{Float64}, d2);

num_trace  = size(d2, 2);
max_delay  = 1.5;
time_delay = max_delay*rand(num_trace);

# apply time shift to the second far shot
dt = 0.004;
d_blend = blend_forward(d1, d2, time_delay, dt);

# deblend by inversion
(m1, m2) = deblending_by_inversion(d_blend, time_delay, dt; max_iter=130, min_value=1e-5)





"""
   From second dimension determine the number of boat
"""
