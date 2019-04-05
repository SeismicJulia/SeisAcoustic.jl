"""
   initialize the pml damping profile in x- and z-direction
"""
function init_pml(nz::Ti, nx::Ti, h::Tv, apml::Tv, npml::Ti,
         omegac::Complex{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    (dampz, dampzb) = damping_profile(nz, npml, apml, h, omegac)
    (dampx, dampxb) = damping_profile(nx, npml, apml, h, omegac)

    return dampz, dampzb, dampx, dampxb
end

"""
   set up damping profile for PML boundary conditions
"""
function damping_profile(n::Ti, npml::Ti, apml::Tv, h::Tv,
                         omegac::Complex{Tv}) where {Ti<:Int64, Tv<:Float64}

    etype = typeof(h)

    # some constants
    alphad    = 1.0
    pi_over_2 = pi * 0.5

    # values of damp and dampb outside of pml equal to one
    damp  = ones(Complex{etype}, n+2)
    dampb = ones(Complex{etype}, n+2)

    xpml  = etype(npml) * h
    xmax  = etype(n -1) * h

    # both side of damp and one side of dampb
    for i = 1 : npml

        x  = (i-1    ) * h
        xb = (i-1+0.5) * h

        epsilon  = apml         * (1.0 - cos((xpml-x )/xpml * pi_over_2))
        epsilonb = apml         * (1.0 - cos((xpml-xb)/xpml * pi_over_2))
        alpha    = (alphad-1.0) * (1.0 - cos((xpml-x )/xpml * pi_over_2)) + 1.0
        alphab   = (alphad-1.0) * (1.0 - cos((xpml-xb)/xpml * pi_over_2)) + 1.0

        damp[ i+1] = 1.0 / (alpha  + im * epsilon  / omegac)
        dampb[i+1] = 1.0 / (alphab + im * epsilonb / omegac)

        # damp profile is symmetrical
        damp[n-i+2] = damp[i+1]

    end
    damp[  1] = damp[2  ]
    damp[n+2] = damp[n+1]

    # right side of dampb
    for i = 1 : npml + 1

        xb  = xmax + 0.5 * h - (i-1)*h

        epsilonb = apml         * (1.0 - cos((xb-xmax+xpml)/xpml * pi_over_2))
        alphab   = (alphad-1.0) * (1.0 - cos((xb-xmax+xpml)/xpml * pi_over_2)) + 1.0

        dampb[n-i+2] = 1.0/(alphab + im * epsilonb / omegac)
    end
    dampb[1  ] = dampb[2]
    dampb[n+2] = dampb[n+1]

    return damp, dampb

end

"""
  compute complex bulk modulus which support visco-acoustic wave-equation
"""
function vel2mu_cmplx(vel::Matrix{Tv}, q::Matrix{Tv}, rho::Matrix{Tv}) where {Tv<:Float64}

    etype = eltype(vel)
    (nz, nx) = size(vel)

    mu = zeros(Complex{etype}, nz+2, nx+2)

    # central part
    for ix = 1 : nx
        for iz = 1 : nz
            vc = vel[iz,ix] * (1.0 - 0.5*im / q[iz,ix])
            mu[iz+1,ix+1] = rho[iz,ix] * vc * vc
        end
    end

    # left and right
    for iz = 1 : nz
        mu[iz+1,1   ] = mu[iz+1,2   ]
        mu[iz+1,nx+2] = mu[iz+1,nx+1]
    end

    # upper and bottom
    for ix = 1 : nx
        mu[1   ,ix+1] = mu[2   ,ix+1]
        mu[nz+2,ix+1] = mu[nz+1,ix+1]
    end

    # four corners
    mu[1   ,1   ] = mu[2   ,2   ]
    mu[1   ,nx+2] = mu[2   ,nx+1]
    mu[nz+2,1   ] = mu[nz+1,2   ]
    mu[nz+2,nx+2] = mu[nz+1,nx+1]

    return mu

end


"""
   compute buoyancy from density
"""
function rho2bu(rho::Matrix{Tv}) where {Tv<:Float64}

    etype    = eltype(rho)
    (nz, nx) = size(rho)

    bu = zeros(etype, nz+2, nx+2)
    for ix = 1 : nx
        for iz = 1 : nz
            if rho[iz  ,ix  ] != 0.0
               bu[ iz+1,ix+1]  = 1.0 / rho[iz,ix]
            else
               bu[ iz+1,ix+1]  = 1e10
            end
        end
    end

    for iz = 1 : nz
        bu[iz+1,1   ] = bu[iz+1,2   ]
        bu[iz+1,nx+2] = bu[iz+1,nx+1]
    end

    for ix = 1 : nx
        bu[1   ,ix+1] = bu[2   ,ix+1]
        bu[nz+2,ix+1] = bu[nz+1,ix+1]
    end

    # four corners
    bu[1   ,1   ] = bu[2   ,2   ]
    bu[1   ,nx+2] = bu[2   ,nx+1]
    bu[nz+2,1   ] = bu[nz+1,2   ]
    bu[nz+2,nx+2] = bu[nz+1,nx+1]

    return bu

end

"""
   get the average of buoyancy
"""
function bu_average(bu::AbstractArray{Tv,2}, i1::Ti, i2::Ti) where {Tv<:Float64, Ti<:Int64}

    i = i1 + 1
    j = i2 + 1

    bmm = 0.25 * (bu[i-1,j-1] + bu[i  ,j-1] + bu[i-1,j  ] + bu[i  ,j  ])
    bpm = 0.25 * (bu[i  ,j-1] + bu[i+1,j-1] + bu[i  ,j  ] + bu[i+1,j  ])

    bmp = 0.25 * (bu[i-1,j  ] + bu[i  ,j  ] + bu[i-1,j+1] + bu[i  ,j+1])
    bpp = 0.25 * (bu[i  ,j  ] + bu[i+1,j  ] + bu[i  ,j+1] + bu[i+1,j+1])

    bm0 = 0.5  * (bu[i  ,j  ] + bu[i-1,j  ]                            )
    bp0 = 0.5  * (bu[i  ,j  ] + bu[i+1,j  ]                            )

    b0m = 0.5  * (bu[i  ,j  ] + bu[i  ,j-1]                            )
    b0p = 0.5  * (bu[i  ,j  ] + bu[i  ,j+1]                            )

    b00 = bu[i,j]

    return bmm, b0m, bpm, bm0, b00, bp0, bmp, b0p, bpp
end

"""
   the kernal of computing the nonzero elements of helmholtz matrix, output is the
   LU factorization of the original matrix.
"""
function helmholtz_kernal(mu::AbstractArray{Complex{Tv}, 2}, bu::AbstractArray{Tv,2}, nz::Ti, nx::Ti,
                          npml::Ti, apml::Tv, h::Tv, omegac::Complex{Tv}) where {Tv<:Float64, Ti<:Int64}

    etype = eltype(bu)

    (dampz, dampzb, dampx, dampxb) = init_pml(nz, nx, h, apml, npml, omegac)

    # some constants for mixed grid combination and lump mass average
    wm1  = 0.6287326      #(c)
    wm2  = 0.3712667      #(4d)
    wm3  = 1 - wm1 - wm2  #(1-c-4d)

    wm2  = 0.25 * wm2     # d
    wm3  = 0.25 * wm3     # 1/4 * (1-c-4d)

    w1   = 0.4382634      # 1-a
    w2   = 1 - w1         # a

    omega2      = omegac * omegac
    one_over_h2 = 1.0 / (h * h)

    nnz = 0
    irn = zeros(Int64         , 9 * nz * nx)
    icn = zeros(Int64         , 9 * nz * nx)
    mat = zeros(Complex{etype}, 9 * nz * nx)

    for i2 = 1 : nx
        for i1 = 1 : nz
            (bmm, b0m, bpm, bm0, b00, bp0, bmp, b0p, bpp) = bu_average(bu, i1, i2)

            i = i1 + 1
            j = i2 + 1

            d2p = dampxb[j  ]
            d2m = dampxb[j-1]
            d1p = dampzb[i  ]
            d1m = dampzb[i-1]

            # node (0, 0)
            row_idx = (i2-1) * nz + i1
            col_idx = row_idx
            nnz = nnz + 1

            irn[nnz] = row_idx
            icn[nnz] = col_idx

            mat[nnz] = ( wm1  * omega2 / mu[i,j]
                        -0.25 * w1 * one_over_h2 * ( dampx[j] * (bmp*d2p + bpm*d2m + bpp*d2p + bmm*d2m)
                                                    +dampz[i] * (bmp*d1m + bpm*d1p + bpp*d1p + bmm*d1m))
                        -1.0  * w2 * one_over_h2 * ( dampx[j] * (b0p*d2p + b0m*d2m)
                                                    +dampz[i] * (bp0*d1p + bm0*d1m)) )

            # node (0, 1)
            if (i2 < nx)
               col_idx = (i2-1+1) * nz + i1
               nnz = nnz + 1
               irn[nnz] = row_idx
               icn[nnz] = col_idx

               mat[nnz] = ( wm2  * omega2 / mu[i,j+1]
                           +0.25 * w1 * one_over_h2 * ( dampx[j] * ( bmp*d2p + bpp*d2p)
                                                       +dampz[i] * (-bmp*d1m - bpp*d1p))
                           +1.0  * w2 * one_over_h2 *   dampx[j] *   b0p*d2p             )
            end

            # node(0, -1)
            if (i2 > 1)
               col_idx = (i2-1-1) * nz + i1
               nnz = nnz + 1
               irn[nnz] = row_idx
               icn[nnz] = col_idx

               mat[nnz] = ( wm2  * omega2 / mu[i,j-1]
                           +0.25 * w1 * one_over_h2 * ( dampx[j] * ( bpm*d2m + bmm*d2m)
                                                       +dampz[i] * (-bpm*d1p - bmm*d1m))
                           +1.0  * w2 * one_over_h2 *   dampx[j] *   b0m*d2m             )
            end

            # node(-1, 0)
            if (i1 > 1)
               col_idx = (i2-1) * nz + i1-1
               nnz = nnz + 1
               irn[nnz] = row_idx
               icn[nnz] = col_idx

               mat[nnz] = ( wm2  * omega2 / mu[i-1,j]
                           +0.25 * w1 * one_over_h2 * ( dampx[j] * (-bmp*d2p - bmm*d2m)
                                                       +dampz[i] * ( bmp*d1m + bmm*d1m))
                           +1.0  * w2 * one_over_h2 *   dampz[i] *   bm0*d1m             )
            end

            # node(1, 0)
            if (i1 < nz)
               col_idx = (i2-1) * nz + i1+1
               nnz = nnz + 1
               irn[nnz] = row_idx
               icn[nnz] = col_idx

               mat[nnz] = ( wm2  * omega2 / mu[i+1,j]
                           +0.25 * w1 * one_over_h2 * ( dampx[j] * (-bpm*d2m - bpp*d2p)
                                                       +dampz[i] * ( bpm*d1p + bpp*d1p))
                           +1.0  * w2 * one_over_h2 *   dampz[i] *   bp0*d1p             )
            end

            # node(-1, 1)
            if (i1 > 1) && (i2 < nx)
               col_idx = (i2-1+1) * nz + i1-1
               nnz = nnz + 1
               irn[nnz] = row_idx
               icn[nnz] = col_idx

               mat[nnz] = ( wm3  * omega2 / mu[i-1,j+1]
                           +0.25 * w1 * one_over_h2 * (dampx[j] * bmp*d2p + dampz[i] * bmp*d1m))
            end

            # node(1,-1)
            if (i1 < nz) && (i2 > 1)
               col_idx = (i2-1-1) * nz + i1+1
               nnz = nnz + 1
               irn[nnz] = row_idx
               icn[nnz] = col_idx

               mat[nnz] = ( wm3  * omega2 / mu[i+1,j-1]
                           +0.25 * w1 * one_over_h2 * (dampx[j] * bpm*d2m + dampz[i] * bpm*d1p))
            end

           # node(1, 1)
           if (i1 < nz) && (i2 < nx)
              col_idx = (i2-1+1) * nz + i1+1
              nnz = nnz + 1
              irn[nnz] = row_idx
              icn[nnz] = col_idx

              mat[nnz] = ( wm3  * omega2 / mu[i+1,j+1]
                          +0.25 * w1 * one_over_h2 * (dampx[j] * bpp*d2p + dampz[i] * bpp*d1p))
           end

           # node(-1,-1)
           if (i1 > 1) && (i2 > 1)
              col_idx = (i2-1-1) * nz + i1-1
              nnz = nnz + 1
              irn[nnz] = row_idx
              icn[nnz] = col_idx

              mat[nnz] = ( wm3  * omega2 / mu[i-1,j-1]
                          +0.25 * w1 * one_over_h2 * (dampx[j] * bmm*d2m + dampz[i] * bmm*d1m))
           end

        end
    end

    H = sparse(irn[1:nnz], icn[1:nnz], mat[1:nnz])
    H = lufact(H)

    return H

end

"""
   For sparse matrix in coordinate form, remove the repeated ones
"""
function remove_repeated_samples(m::Ti, n::Ti,
         irow::Ti, icol::Ti, nzval::Tv) where {Ti<:Int64, Tv<:Number}

  etype = eltype(nzval)
  nnz   = length(nzval)

  colptr = zeros(Int64, n+1)
  rowidx = zeros(Int64, nnz)
  nzval1 = zeros(etype, nnz)

  # Determine the column length
  for i = 1 : nnz
      colptr[icol[i]] = colptr[icol[i]] + 1
  end

  # starting position of each column
  k = 1
  for j = 1 : N + 1
      k0 = colptr[j]
      colptr[j] = k
      k = k + k0
  end

  # fill in the output
  for k = 1 : nnz
      i = irow[k]
      j = icol[k]
      x = nzval[k]
      idx = colptr[j]
      nzval1[idx] = x
      rowidx[idx] = i
      colptr[j  ] = idx + 1
  end

  # shift back colptr
  for j = N : -1 : 1
      colptr[j+1] = colptr[j]
  end
  colptr[1] = 1

  A = spzeros(m, n)
  A.rowval = rowidx
  A.colptr = colptr
  A.nzval  = nzval1

  return A

end

"""
   convert sparse matrix from coordinate format to CSC format
"""

function COOCSC(m::Ti, n::Ti,
         irow::Vector{Ti}, icol::Vector{Ti}, nzval::Vector{Tv}) where {Ti<:Int64, Tv<:Number}

  etype = eltype(nzval)
  nnz   = length(nzval)

  colptr = zeros(Int64, n+1)
  rowidx = zeros(Int64, nnz)
  nzval1 = zeros(etype, nnz)

  # Determine the column length
  for i = 1 : nnz
      colptr[icol[i]] = colptr[icol[i]] + 1
  end

  # starting position of the element in each column
  k = 1
  for j = 1 : n + 1
      k0 = colptr[j]
      colptr[j] = k
      k = k + k0
  end

  # fill in the output
  for k = 1 : nnz
      i = irow[k]
      j = icol[k]
      x = nzval[k]
      idx = colptr[j]
      nzval1[idx] = x
      rowidx[idx] = i
      colptr[j  ] = idx + 1
  end

  # shift back colptr
  for j = n : -1 : 1
      colptr[j+1] = colptr[j]
  end
  colptr[1] = 1

  A = spzeros(m, n)
  A.rowval = rowidx
  A.colptr = colptr
  A.nzval  = nzval1

  return A

end

"""
   evaluate the helmholtz matrix for one frequency and do LU decomposition.
"""
function eval_helmholtz_matrix(vel::Matrix{Tv}, q::Matrix{Tv}, rho::Matrix{Tv},
                               nzm::Ti, nxm::Ti, npml::Ti, apml::Tv, h::Tv,
                               omegac::Complex{Tv}) where {Ti<:Int64, Tv<:Float64}

    etype = eltype(vel)
    nz = nzm + 2 * npml
    nx = nxm + 2 * npml

    # model parameter extending
    vele = model_extend(vel, npml)
    qe   = model_extend(q  , npml)
    rhoe = model_extend(rho, npml)
    mu   = vel2mu_cmplx(vele, qe, rhoe)
    bu   = rho2bu(rhoe)

    H    = helmholtz_kernal(mu, bu, nz, nx, npml, apml, h, omegac)

    return H

end


"""
   Generate wavefield in time domain by solving acoustic wave equation in frequency domain
"""
function eval_wavefield_FDFD(nzm::Ti, nxm::Ti, npml::Ti, apml::Tv, h::Tv,
                             vel::Matrix{Tv}, q::Matrix{Tv}, rho::Matrix{Tv},
                             isz::Ti, isx::Ti, wavelet::Vector{Tv},
                             dt::Tv, tmax::Tv, flower::Tv, fupper::Tv; print_flag=false) where {Ti<:Int64, Tv<:Float64}

    etype = eltype(vel)

    nz = nzm + 2*npml
    nx = nxm + 2*npml

    vele = model_extend(vel, npml)
    qe   = model_extend(q  , npml)
    rhoe = model_extend(rho, npml)
    mu   = vel2mu_cmplx(vele, qe, rhoe)
    bu   = rho2bu(rhoe)

    # make the number of time samples to be even
    # good to take advantage complex conjugate property
    nt = round(Int64, tmax/dt) + 1
    if mod(nt, 2) != 0
       nt = nt + 1
    end

    # allocate space for wavefield
    # first dimension is time
    wcube = zeros(Complex128, nt, nz, nx)

    # maximum of frequency
    fmax = 1.0  / dt
    df   = fmax / nt

    iw_lower = floor(Int64, flower / df) + 1
    iw_upper = ceil( Int64, fupper / df) + 1
    if print_flag
       println("iw_lower=$iw_lower")
       println("iw_upper=$iw_upper")
    end

    # transform source wavelet to frequency domain
    if length(wavelet) < nt
       wavelet = vcat(wavelet, zeros(nt-length(wavelet)))
    elseif length(wavelet) > nt
       error("the source wavelet are too long")
    end

    # opposite convention of fourier transform
    fwave = sqrt(nt) * ifft(wavelet)

    # one frequency sample of source
    b = zeros(Complex128, nz, nx)
    isz_pad = isz + npml
    isx_pad = isx + npml

    # loop over frequency slice
    for iw = iw_lower : iw_upper

        omegac = Complex128(2.0 * pi * (iw-1) * df)
        H = helmholtz_kernal(mu, bu, nz, nx, npml, apml, h, omegac)

        #the correct way to implement source
        b[isz_pad, isx_pad] = im * omegac * fwave[iw] / (vel[isz,isx]*vel[isz,isx])
        x  = H \ vec(b)

        # symmetry property
        tmp= reshape(x, nz, nx)
        wcube[iw     ,:,:] = tmp
        wcube[nt-iw+2,:,:] = conj(tmp)

        if print_flag
           println("finished $iw")
        end
    end

    # transform back to time domain
    wcube = real(1/sqrt(nt) * fft(wcube, 1))

    return wcube

end
