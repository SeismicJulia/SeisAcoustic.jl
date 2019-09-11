function laplace_filter(x::Matrix{Tv}) where {Tv<:AbstractFloat}

    (n1, n2) = size(x)
    x1 = zeros(Tv, n1, n2)
    for i2 = 2 : n2-1
        for i1 = 2 : n1-1
            x1[i1,i2] = (8.0 * x[i1,i2] - x[i1-1,i2-1] - x[i1,i2-1] - x[i1+1,i2-1]
                                        - x[i1-1,i2  ]              - x[i1+1,i2  ]
                                        - x[i1-1,i2+1] - x[i1,i2+1] - x[i1+1,i2+1])
        end
    end
    return x1
end

# function laplace_filter(img::Matrix{Tv}, iflag::Int64) where {Tv<:AbstractFloat}
#     (m, n) = size(img)
#     # d = zeros(eltype(img), m, n)
#     # for ix = 2 : n-1
#     #     for iz = 2 : m-1
#     #         d[iz, ix] = 4*img[iz,ix]-img[iz-1,ix]-img[iz+1,ix]-img[iz,ix-1]-img[iz,ix+1]
#     #     end
#     # end
#     a = -1.0 * ones(m-1); a[end] = 0.0
#     b =    2 * ones(m  ); b[1]   = 1.0; b[end] = 1.0
#     c = -1.0 * ones(m-1); c[1]   = 0.0;
#     V = spdiagm((a,b,c), (-1,0,1), m, m)
#
#     a = -1.0 * ones(n-1); a[end] = 0.0
#     b =    2 * ones(n  ); b[1]   = 1.0; b[end] = 1.0
#     c = -1.0 * ones(n-1); c[1]   = 0.0;
#     H = spdiagm((a,b,c), (-1,0,1), n, n)
#
#     if iflag == 1
#        r = V * img + img * H'
#     elseif iflag == 2
#        r = V'* img + img * H
#     end
#
#     return r
# end
