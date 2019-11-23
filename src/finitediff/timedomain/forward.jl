"""
   compute the partial derivative in vertical direction
"""
function vertical_partial_derivative!(r::Vector{Tv}, p::Vector{Tv}, dpdz::SparseMatrixCSC{Tv,Ti},
         nz::Ti, nx::Ti, tmp1::Vector{Tv}, tmp2::Vector{Tv}; iflag=1) where {Ti<:Int64, Tv<:AbstractFloat}

    istart = 1
    for ix = 1 : nx

        # process column by column
        copyto!(tmp1, 1, p, istart, nz)

        # main part
        if iflag == 1
           A_mul_b!(tmp2, dpdz, tmp1)

        elseif iflag == 2
           At_mul_b!(tmp2, dpdz, tmp1)
        end

        # save the result
        copyto!(r, istart, tmp2, 1, nz)

        # prepare for next column
        istart = istart + nz
    end

    return nothing
end

"""
   compute the partial derivative in horizontal direction
"""
function horizontal_partial_derivative!(r::Vector{Tv}, p::Vector{Tv}, dpdx::SparseMatrixCSC{Tv,Ti},
         nz::Ti, nx::Ti, tmp1::Vector{Tv}, tmp2::Vector{Tv}; iflag=1) where {Ti<:Int64, Tv<:AbstractFloat}

    # the address of vector p and r
    sp = pointer(p)
    sr = pointer(r)

    # electronic size
    esize = sizeof(Tv)
    shift = 0

    for iz = 1 : nz

        # extract one row
        BLAS.blascopy!(nx, sp+shift, nz, tmp1, 1)

        # main part
        if iflag == 1
           A_mul_b!(tmp2, dpdx, tmp1)

        elseif iflag == 2
           At_mul_b!(tmp2, dpdx, tmp1)
        end

        # save the result
        BLAS.blascopy!(nx, tmp2, 1, sr+shift, nz)

        # prepare for next row
        shift = shift + esize
    end

    return nothing
end

"""
   one step forward modelling
"""
function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, params::TdParams{Ti,Tv},
         tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
         tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    for ix = 1 : params.Nx
        ilower = (ix-1) * params.Nz + 1
        iupper = ilower + params.Nz - 1

        # copy one column
        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = spt1.pz[iz] + spt1.px[iz]  # p1 = pz1 + px1
        end

        # partial derivative in Z direction
        A_mul_b!(tmp_z2, params.dpdz, tmp_z1)        # dpdz
        # mul!(tmp_z2, params.dpdz, tmp_z1)

        # update one column
        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            spt2.vz[iz] = params.MvzBvz[idx]*spt1.vz[iz] + params.MvzBp[iz]*tmp_z2[idx]
        end
    end

    # update vx row-by-row
    for iz = 1 : params.Nz

        idx = iz
        for ix = 1 : params.Nx
            tmp_x1[ix] = spt1.pz[idx] + spt1.px[idx] # p1 = pz1 + px1
            idx= idx + params.Nz
        end

        A_mul_b!(tmp_x2, params.dpdx, tmp_x1)        # dpdx
        # mul!(tmp_x2, params.dpdx, tmp_x1)        # dpdx

        idx = iz
        for ix = 1 : params.Nx
            spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*tmp_x2[ix]
            idx= idx + params.Nz
        end
    end

    # update pz column-by-column
    for ix = 1 : params.Nx
        ilower = (ix-1) * params.Nz + 1
        iupper = ilower + params.Nz - 1

        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = spt2.vz[iz]
        end

        A_mul_b!(tmp_z2, params.dvdz, tmp_z1)        # dpdz
        # mul!(tmp_z2, params.dvdz, tmp_z1)        # dpdz


        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            spt2.pz[iz] = params.MpzBpz[idx]*spt1.pz[iz] + params.MpzBvz[iz]*tmp_z2[idx]
        end
    end

    # update px row-by-row
    for iz = 1 : params.Nz

        idx = iz
        for ix = 1 : params.Nx
            tmp_x1[ix] = spt2.vx[idx]
            idx= idx + params.Nz
        end

        A_mul_b!(tmp_x2, params.dvdx, tmp_x1)         # dvdx
        # mul!(tmp_x2, params.dvdx, tmp_x1)         # dvdx

        idx = iz
        for ix = 1 : params.Nx
            spt2.px[idx] = params.MpxBpx[ix]*spt1.px[idx] + params.MpxBvx[idx]*tmp_x2[ix]
            idx= idx + params.Nz
        end
    end

    return nothing
end

# ==============================================================================
#                        forward simulation
# ==============================================================================
"""
   forward modelling and sampling the wavefield at locations specified by (rz, rx),
the shot gather can also be written to "path_shot"
"""
function multi_step_forward!(rz::Vector, rx::Vector, src::Ts, params::TdParams;
         location_flag="index", path_shot="NULL") where {Ts<:Union{Source, Vector{Source}}}

    # allocate space for recordings
    rec = Recordings(rz, rx, params; location_flag=location_flag)

    # initialize some variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)
    sample_spt2rec!(rec, spt1, 1)

    # applying time stepping nt-1 times
    for it = 2 : params.nt

        # one step forward
        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        add_source!(spt2, src, it)
        sample_spt2rec!(rec, spt2, it)

        # prepare for next time stepping
        copy_snapshot!(spt1, spt2)
    end

    if path_shot == "NULL"
       return rec
       # return Recordings(rec.nt, rec.nr, rec.dt, rec.irz, rec.irx, rec.spt2rec,
       #        convert(Matrix{params.data_format}, rec.p))
    else
       write_recordings(path_shot, rec)
       return path_shot
    end
end

"""
   save snapshot, full wave field and pressure field, which are generated by a single source, to disk
path_spt != "NULL" : vz, vx, pz, px include PML part, every save_interval step
path_wfd != "NULL" : vz, vx, pz+px  without PML part, every save_interval step
path_pre != "NULL" : p=pz+px without PML bounary part. every save_interval step
path_bnd != "NULL" : save boundary of wavefield at every step
path_lwfd!= "NULL" : save the last wavefield
path_sws != "NULL" : save the source-side wavefield strength (preconditioner for RTM)
"""
function multi_step_forward!(src::Ts, params::TdParams;
                             path_spt="NULL", path_wfd="NULL" , path_pre="NULL", interval=200,
                             path_bnd="NULL", path_lwfd="NULL", path_sws="NULL") where {Ts<:Union{Source, Vector{Source}}}

    # initialize some variables
    spt1   = Snapshot(params)
    spt2   = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # save snapshot
    if path_spt != "NULL"
       hdr_spt   = snapshot_header(params, interval)
       fid_spt   = write_RSheader(path_spt, hdr_spt)
       append_one_snapshot(fid_spt, spt1)
    end

    # save wavefield
    if path_wfd != "NULL"
       hdr_wfd   = wavefield_header(params, interval)
       fid_wfd   = write_RSheader(path_wfd, hdr_wfd)
       append_one_wavefield(fid_wfd, spt1, params)
    end

    # save pressure
    if path_pre != "NULL"
       hdr_pre   = pressure_header(params, interval)
       fid_pre   = write_RSheader(path_pre, hdr_pre)
       append_one_pressure(fid_pre, spt1, params)
    end

    # boundary of wavefield
    if path_bnd != "NULL"
       hdr_bnd = boundary_header(params)
       fid_bnd = write_RSheader(path_bnd, hdr_bnd)
       append_one_boundary(fid_bnd, spt1, params)
    end

    # strength of source-side wavefield
    if path_sws != "NULL"
       N = params.nz * params.nx
       strength = zeros(params.data_format, N)
    end

    # applying time stepping nt-1 times
    for it = 2 : params.nt

        # one step forward
        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # compute source-side wavefield
        path_sws != "NULL" && source_strength!(strength, spt1, spt2, params)

        add_source!(spt2, src, it)

        # write current wavefield to disk
        if mod(it-1, interval) == 0
           path_spt != "NULL" && append_one_snapshot(fid_spt, spt2)
           path_wfd != "NULL" && append_one_wavefield(fid_wfd, spt2, params)
           path_pre != "NULL" && append_one_pressure(fid_pre, spt2, params)
        end

        # save boundary of wavefield
        path_bnd != "NULL" && append_one_boundary(fid_bnd, spt2, params)

        # prepare for next time stepping
        copy_snapshot!(spt1, spt2)
    end

    # save the source-side wavefield strength to disk
    if path_sws != "NULL"
       strength  = reshape(strength, params.nz, params.nx)
       hdr_sws   = RegularSampleHeader(strength; title="source strength")
       write_RSdata(path_sws, hdr_sws, strength)
    end

    # save the last wavefield
    if path_lwfd != "NULL"
       wfd = sample_spt2wfd(spt1, params)
       write_wavefield(path_lwfd, wfd, params)
    end

    # close file
    path_spt != "NULL" && close(fid_spt)
    path_wfd != "NULL" && close(fid_wfd)
    path_pre != "NULL" && close(fid_pre)
    path_bnd != "NULL" && close(fid_bnd)

    return nothing
end

# define a auxillary function to improve efficiency
function source_strength!(strength, spt1, spt2, params)

    # total number of elements
    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]

        p2= spt2.pz[j]+spt2.px[j]
        p1= spt1.pz[j]+spt1.px[j]

        tmp = 2.0 * (p2 - p1) / params.vel[i]
        strength[i] += tmp * tmp
    end
    return nothing
end

# ==============================================================================
#          one step forward modelling with saved boundary value
# ==============================================================================
"""
   one step-forward reconstruction of wavefield via recording the boundary value
"""
function one_step_forward!(wfd2::Wavefield, wfd1::Wavefield,
         bnd::WavefieldBound, params::TdParams,
         tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
         tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    # update vz column-by-column
    for ix = 1 : params.nx
        ilower = (ix-1) * params.nz + 1
        iupper = ilower + params.nz - 1

        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = wfd1.p[iz]
        end

        A_mul_b!(tmp_z2, params.rpdz, tmp_z1)        # dpdz

        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            wfd2.vz[iz] = wfd1.vz[iz] + params.RvzBp[iz]*tmp_z2[idx]
        end
    end

    # update vx row-by-row
    for iz = 1 : params.nz

        idx = iz
        for ix = 1 : params.nx
            tmp_x1[ix] = wfd1.p[idx]
            idx= idx + params.nz
        end

        A_mul_b!(tmp_x2, params.rpdx, tmp_x1)         # dpdx

        idx = iz
        for ix = 1 : params.nx
            wfd2.vx[idx] = wfd1.vx[idx] + params.RvxBp[idx]*tmp_x2[ix]
            idx= idx + params.nz
        end
    end

    # correct for boundary part
    for i = 1 : length(params.wfd2bnd)
        j = params.wfd2bnd[i]
        wfd2.vz[j] = bnd.vz[i]
        wfd2.vx[j] = bnd.vx[i]
    end

    # update pz column-by-column
    for ix = 1 : params.nx
        ilower = (ix-1) * params.nz + 1
        iupper = ilower + params.nz - 1

        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = wfd2.vz[iz]
        end

        A_mul_b!(tmp_z2, params.rvdz, tmp_z1)        # dpdz

        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            wfd2.p[iz] = wfd1.p[iz] + params.RpBv[iz]*tmp_z2[idx]
        end
    end

    # update px row-by-row
    for iz = 1 : params.nz

        idx = iz
        for ix = 1 : params.nx
            tmp_x1[ix] = wfd2.vx[idx]
            idx= idx + params.nz
        end

        A_mul_b!(tmp_x2, params.rvdx, tmp_x1)         # dvdx

        idx = iz
        for ix = 1 : params.nx
            wfd2.p[idx] = wfd2.p[idx] + params.RpBv[idx]*tmp_x2[ix]
            idx= idx + params.nz
        end
    end

    # correct for boundary part
    for i = 1 : length(params.wfd2bnd)
        j = params.wfd2bnd[i]
        wfd2.p[j] = bnd.p[i]
    end

    return nothing
end

"""
   Forward reconstruct acoustic pressure field using the saved boundary wavefield value and
the source (used for concept proofing)
"""
function pressure_reconstruct_forward(path_bnd::Tp, src::Ts,
         params::TdParams) where {Tp<:String, Ts<:Union{Source,Vector{Source}}}

    # make the source as a vector has one element
    if Ts <: Source
       srcs    = Vector{Source}(undef,1)
       srcs[1] = src
    else
       srcs    = src
    end
    ns = length(srcs)

    # length of one-step pressure field
    N = params.nz * params.nx

    # initialize intermediate variables
    wfd1   = Wavefield(params)
    wfd2   = Wavefield(params)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd = WavefieldBound(params)

    # allocate memory for saving pressure field
    pre = zeros(params.data_format, N * params.nt)

    # add source to the first wavefield
    add_source!(wfd1, srcs, 1)

    # save pressure field
    copyto!(pre, 1, wfd1.p, 1, N)

    # the start index for saving the next pressure field
    idx_o = N+1

    # loop over time stepping
    for it = 2 : params.nt

        # read the boundary value
        read_one_boundary!(bnd, fid_bnd, it, params)

        # forward time steping and correcting boundaries
        one_step_forward!(wfd2, wfd1, bnd, params,
                          tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # loop over sources
        for i = 1 : ns
            if (params.order+1 <= srcs[i].isz <= params.nz-params.order &&
                params.order+1 <= srcs[i].isx <= params.nx-params.order)

               add_source!(wfd2, srcs[i], it)
            end
        end

        # save the pressure field
        copyto!(pre, idx_o, wfd2.p, 1, N)
        idx_o = idx_o + N

        # prepare for next step
        copy_wavefield!(wfd1, wfd2)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt)
end

# this version of one step forward is little bit slower than the used one
# since the overhead of calling an function
# function one_step_forward_new!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, params::TdParams{Ti,Tv},
#                                tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
#                                tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # compute pressure field
#     spt2.pz .= spt1.pz .+ spt1.px
#
#     # update particle velocity
#     vertical_partial_derivative!(spt2.vz, spt2.pz, params.dpdz, params.Nz, params.Nx, tmp_z1, tmp_z2; iflag=1)
#     horizontal_partial_derivative!(spt2.vx, spt2.pz, params.dpdx, params.Nz, params.Nx, tmp_x1, tmp_x2; iflag=1)
#     idx = 0
#     for ix = 1 : params.Nx
#         for iz = 1 : params.Nz
#             idx= idx + 1
#             spt2.vz[idx] = params.MvzBvz[iz]*spt1.vz[idx] + params.MvzBp[idx]*spt2.vz[idx]
#             spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*spt2.vx[idx]
#         end
#     end
#
#     # update pressure field
#     vertical_partial_derivative!(spt2.pz, spt2.vz, params.dvdz, params.Nz, params.Nx, tmp_z1, tmp_z2; iflag=1)
#     horizontal_partial_derivative!(spt2.px, spt2.vx, params.dvdx, params.Nz, params.Nx, tmp_x1, tmp_x2; iflag=1)
#     idx = 0
#     for ix = 1 : params.Nx
#         for iz = 1 : params.Nz
#             idx= idx + 1
#             spt2.pz[idx] = params.MpzBpz[iz]*spt1.pz[idx] + params.MpzBvz[idx]*spt2.pz[idx]
#             spt2.px[idx] = params.MpxBpx[ix]*spt1.px[idx] + params.MpxBvx[idx]*spt2.px[idx]
#         end
#     end
#
#     return nothing
# end

# function one_step_forward!(wfd2::Wavefield, wfd1::Wavefield,
#          bnd::WavefieldBound, params::TdParams,
#          tmp_a1::Vector{Tv}, tmp_a2::Vector{Tv},
#          tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
#          tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # total number of field elements
#     N = params.nz * params.nx
#
#     # partial derivative in vertical direction
#     vertical_partial_derivative!(tmp_a1, wfd1.p, params.rpdz,
#     params.nz, params.nx, tmp_z1, tmp_z2)
#
#     # partial derivative in horizontal direction
#     horizontal_partial_derivative!(tmp_a2, wfd1.p, params.rpdx,
#     params.nz, params.nx, tmp_x1, tmp_x2)
#
#     # update vz
#     for i = 1 : N
#         wfd2.vz[i] = wfd1.vz[i] + params.RvzBp[i]*tmp_a1[i]
#         wfd2.vx[i] = wfd1.vx[i] + params.RvxBp[i]*tmp_a2[i]
#     end
#
#     # correct for boundary part
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd2.vz[j] = bnd.vz[i]
#         wfd2.vx[j] = bnd.vx[i]
#     end
#
#     # partial derivative in vertical direction
#     vertical_partial_derivative!(tmp_a1, wfd2.vz, params.rvdz,
#     params.nz, params.nx, tmp_z1, tmp_z2)
#
#     # partial derivative in horizontal direction
#     horizontal_partial_derivative!(tmp_a2, wfd2.vx, params.rvdx,
#     params.nz, params.nx, tmp_x1, tmp_x2)
#
#     # update for pressure
#     for i = 1 : N
#         wfd2.p[i] = wfd1.p[i] + params.RpBv[i] * (tmp_a1[i] + tmp_a2[i])
#     end
#
#     # correct for boundary part
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd2.p[j] = bnd.p[i]
#     end
#
#     return nothing
# end
