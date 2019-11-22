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
    wfd1 = Wavefield(params)
    wfd2 = Wavefield(params)
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


# function multi_step_forward!(rec::Recordings, srcs::Vector{Source}, params::TdParams;
#                              path_spt="NULL", path_wfd="NULL", path_pre="NULL", interval=1)
#
#     # initialize some variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp_z1 = zeros(params.data_format, params.Nz)
#     tmp_z2 = zeros(params.data_format, params.Nz)
#     tmp_x1 = zeros(params.data_format, params.Nx)
#     tmp_x2 = zeros(params.data_format, params.Nx)
#
#     # add source to the first snapshot
#     add_multi_sources!(spt1, srcs, 1)
#     sample_spt2rec!(rec, spt1, 1)
#
#     # save snapshot
#     if path_spt != "NULL"
#        hdr_spt   = snapshot_header(params, interval)
#        fid_spt   = write_RSheader(path_spt, hdr_spt)
#        append_one_snapshot(fid_spt, spt1)
#     end
#
#     # save wavefield
#     if path_wfd != "NULL"
#        hdr_wfd   = wavefield_header(params, interval)
#        fid_wfd   = write_RSheader(path_wfd, hdr_wfd)
#        append_one_wavefield(fid_wfd, spt1, params)
#     end
#
#     # save pressure
#     if path_pre != "NULL"
#        hdr_pre   = pressure_header(params, interval)
#        fid_pre   = write_RSheader(path_pre, hdr_pre)
#        append_one_pressure(fid_pre, spt1, params)
#     end
#
#     # applying time stepping nt-1 times
#     for it = 2 : params.nt
#
#         # one step forward
#         one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
#         add_multi_sources!(spt2, srcs, it)
#         sample_spt2rec!(rec, spt2, it)
#
#         # write current wavefield to disk
#         if mod(it-1, interval) == 0
#            path_spt != "NULL" && append_one_snapshot(fid_spt, spt2)
#            path_wfd != "NULL" && append_one_wavefield(fid_wfd, spt2, params)
#            path_pre != "NULL" && append_one_pressure(fid_pre, spt2, params)
#         end
#
#         # prepare for next time stepping
#         copy_snapshot!(spt1, spt2)
#     end
#
#     # close file
#     path_spt != "NULL" && close(fid_spt)
#     path_wfd != "NULL" && close(fid_wfd)
#     path_pre != "NULL" && close(fid_pre)
#
#     return nothing
# end

# """
#    compute the recordings via forward modelling with simultaneouse source, the boundary and the wavefield
# at the last time step are saved if "path_bnd" and "path_wfd" are given.
# """
# function multi_step_forward!(srcs::Vector{Source}, params::TdParams;
#                              path_bnd="NULL", path_wfd="NULL", path_sws="NULL")
#
#     # at least one output
#     path_bnd == "NULL" && path_bnd == "NULL" && path_bnd == "NULL" && error("at least one output")
#
#     # initialize variables for time stepping
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp_z1 = zeros(params.data_format, params.Nz)
#     tmp_z2 = zeros(params.data_format, params.Nz)
#     tmp_x1 = zeros(params.data_format, params.Nx)
#     tmp_x2 = zeros(params.data_format, params.Nx)
#
#     # add source to the first snapshot
#     add_multi_sources!(spt1, srcs, 1)
#
#     # save the boundary of the first wavefield
#     if path_bnd != "NULL"
#        hdr = boundary_header(params)
#        fid = write_RSheader(path_bnd, hdr)
#        append_one_boundary(fid, spt1, params)
#     end
#
#     # save the strength of source-side wavefield
#     if path_sws != "NULL"
#        N = params.nz * params.nx
#        strength = zeros(params.data_format, N)
#     end
#
#     # loop over time stepping
#     for it = 2 : params.nt
#
#         # one time stepping
#         one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
#
#         # compute source-side wavefield
#         path_sws != "NULL" && source_strength!(strength, spt1, spt2, params)
#
#         add_multi_sources!(spt2, srcs, it)
#
#         # save boundary of wavefield
#         path_bnd != "NULL" && append_one_boundary(fid, spt2, params)
#
#         # prepare for next iteration
#         copy_snapshot!(spt1, spt2)
#     end
#
#     # finish saving all the boundary values, close the file
#     path_bnd != "NULL" && close(fid)
#
#     # save the source-side wavefield strength to disk
#     if path_sws != "NULL"
#        strength  = reshape(strength, params.nz, params.nx)
#        hdr       = RegularSampleHeader(strength; title="source strength")
#        write_RSdata(path_sws, hdr, strength)
#     end
#
#     # save the wavefield at the final time step
#     if path_wfd != "NULL"
#        wfd = sample_spt2wfd(spt1, params)
#        write_wavefield(path_wfd, wfd, params)
#     end
#
#     return nothing
# end

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

# function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, params::TdParams{Ti,Tv},
#          tmp_a1::Vector{Tv}, tmp_a2::Vector{Tv},
#          tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
#          tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # sum pz and px
#     copyto!(tmp_a1, spt1.pz)
#     axpy!(1.0, spt1.px, tmp_a1)
#
#     # partial derivative in vertical direction
#     vertical_partial_derivative!(tmp_a2, tmp_a1, params.dpdz,
#     params.Nz, params.Nx, tmp_z1, tmp_z2)
#
#     # update vz
#     idx = 0
#     for ix = 1 : params.Nx
#         for iz = 1 : params.Nz
#             idx = idx + 1
#             spt2.vz[idx] = params.MvzBvz[iz]*spt1.vz[idx] + params.MvzBp[idx]*tmp_a2[idx]
#         end
#     end
#
#     # partial derivative in vertical direction
#     horizontal_partial_derivative!(tmp_a2, tmp_a1, params.dpdx,
#     params.Nz, params.Nx, tmp_x1, tmp_x2)
#
#     # update vx
#     idx = 0
#     for ix = 1 : params.Nx
#         for iz = 1 : params.Nz
#             idx = idx + 1
#             spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*tmp_a2[idx]
#         end
#     end
#
#     # partial derivative of vz in vertical direction
#     vertical_partial_derivative!(tmp_a1, spt2.vz, params.dvdz,
#     params.Nz, params.Nx, tmp_z1, tmp_z2)
#
#     # partial derivative of vx in horizontal direction
#     horizontal_partial_derivative!(tmp_a2, spt2.vx, params.dvdx,
#     params.Nz, params.Nx, tmp_x1, tmp_x2)
#
#     # update pz and px
#     idx = 0
#     for ix = 1 : params.Nx
#         for iz = 1 : params.Nz
#             idx = idx + 1
#             spt2.pz[idx] = params.MpzBpz[iz]*spt1.pz[idx] + params.MpzBvz[idx]*tmp_a1[idx]
#             spt2.px[idx] = params.MpxBpx[ix]*spt1.px[idx] + params.MpxBvx[idx]*tmp_a2[idx]
#         end
#     end
#
#     return nothing
# end

# """
#    compute recordings with a single source
# """
# function multi_step_forward!(rec::Recordings, src::Source, params::TdParams)
#
#     # initialize intermediate variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp_z1 = zeros(params.data_format, params.Nz)
#     tmp_z2 = zeros(params.data_format, params.Nz)
#     tmp_x1 = zeros(params.data_format, params.Nx)
#     tmp_x2 = zeros(params.data_format, params.Nx)
#
#     # add source to the first snapshot
#     add_source!(spt1, src, 1)
#
#     # obtain the first time sample of recordings
#     sample_spt2rec!(rec, spt1, 1)
#
#     # loop over time stepping
#     for it = 2 : params.nt
#
#         one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
#         add_source!(spt2, src, it)
#         sample_spt2rec!(rec, spt2, it)
#
#         # prepare for next time step
#         copy_snapshot!(spt1, spt2)
#     end
#
#     return nothing
# end
#
# """
#    compute recordings with multiple simultaneous source, the only difference is that the input
# is a vector of source
# """
# function multi_step_forward!(rec::Recordings, srcs::Vector{Source}, params::TdParams)
#
#     # initialize intermediate variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp_z1 = zeros(params.data_format, params.Nz)
#     tmp_z2 = zeros(params.data_format, params.Nz)
#     tmp_x1 = zeros(params.data_format, params.Nx)
#     tmp_x2 = zeros(params.data_format, params.Nx)
#
#     # add multi-sources to the first snapshot
#     add_multi_sources!(spt1, srcs, 1)
#
#     # sample the snapshot to fill recordings
#     sample_spt2rec!(rec, spt1, 1)
#
#     # loop over time
#     for it = 2 : params.nt
#
#         one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
#         add_multi_sources!(spt2, srcs, it)
#         sample_spt2rec!(rec, spt2, it)
#
#         # prepare for the next time step
#         copy_snapshot!(spt1, spt2)
#     end
#
#     return nothing
# end

# function multiplication!(c::Vector{Tv}, a::Vector{Tv}, b::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     for i = 1 : length(a)
#         @inbounds c[i] = a[i] * b[i]
#     end
#     return nothing
# end
#
# function multiplication!(b::Vector{Tv}, a::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     for i = 1 : length(a)
#         @inbounds b[i] = b[i] * a[i]
#     end
#     return nothing
# end
#
# function addition!(c::Vector{Tv}, a::Vector{Tv}, b::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     for i = 1 : length(a)
#         @inbounds c[i] = a[i]+b[i]
#     end
#     return nothing
# end
#
# function addition!(a::Vector{Tv}, b::Vector{Tv}) where {Tv<:AbstractFloat}
#     for i = 1 : length(a)
#         @inbounds a[i] = a[i]+b[i]
#     end
#     return nothing
# end
#
# function minus!(a::Vector{Tv}, b::Vector{Tv}, c::Vector{Tv}) where {Tv<:AbstractFloat}
#     for i = 1 : length(a)
#         @inbounds a[i] = b[i] - c[i]
#     end
#     return nothing
# end
#
# """
#    one forward time-stepping of finite-difference
# (rarely used as it take a large amount of memory allocations)
# """
# function one_step_forward!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil)
#
#     spt2.vz = ofds.MvzBvz .* spt1.vz + ofds.MvzBp  * (spt1.pz+spt1.px)
#     spt2.vx = ofds.MvxBvx .* spt1.vx + ofds.MvxBp  * (spt1.pz+spt1.px)
#     spt2.pz = ofds.MpzBpz .* spt1.pz + ofds.MpzBvz *  spt2.vz
#     spt2.px = ofds.MpxBpx .* spt1.px + ofds.MpxBvx *  spt2.vx
#
#     return nothing
# end
#
# """
#    one forward time-stepping of finite-difference
# """
# function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, ofds::ObsorbFDStencil{Ti,Tv},
#          tmp1::Vector{Tv}, tmp2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     addition!(tmp1, spt1.pz, spt1.px)               # p1 = pz1+px1
#
#     # update vz
#     multiplication!(spt2.vz, ofds.MvzBvz, spt1.vz)  # vz2 = MvzBvz * vz1
#     A_mul_b!(tmp2, ofds.MvzBp, tmp1)                # tmp2= MvzBp  * p1
#     addition!(spt2.vz, tmp2)                        # vz2 = vz2    + tmp2
#
#     # update vx
#     multiplication!(spt2.vx, ofds.MvxBvx, spt1.vx)  # vx2 = MvxBvx * vx1
#     A_mul_b!(tmp2, ofds.MvxBp, tmp1)                # tmp2= MvxBp  * p1
#     addition!(spt2.vx, tmp2)                        # vx2 = vx2    + tmp2
#
#     # update pz
#     multiplication!(spt2.pz, ofds.MpzBpz, spt1.pz)  # pz2 = MpzBpz * pz1
#     A_mul_b!(tmp1, ofds.MpzBvz, spt2.vz)            # tmp1= MpzBvz * vz2
#     addition!(spt2.pz, tmp1)                        # pz2 = pz2    + tmp1
#
#     # update px
#     multiplication!(spt2.px, ofds.MpxBpx, spt1.px)  # px2 = MpxBpx * px1
#     A_mul_b!(tmp1, ofds.MpxBvx, spt2.vx)            # tmp1= MpxBvx * vx2
#     addition!(spt2.px, tmp1)                        # px2 = px2    + tmp1
#
#     return nothing
# end
#
# """
#    one step Forward reconstruction of wave field from the boundary values
# """
# function one_step_forward!(wfd2::Wavefield, wfd1::Wavefield, rfds::RigidFDStencil, it::Ti,
#          bnd::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, params::TdParams) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # update vz and vx
#     A_mul_b!(tmp1, rfds.MvzBp, wfd1.p)
#     A_mul_b!(tmp2, rfds.MvxBp, wfd1.p)
#     addition!(wfd2.vz, wfd1.vz, tmp1)
#     addition!(wfd2.vx, wfd1.vx, tmp2)
#
#     # correct for boundary part
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd2.vz[j] = bnd.vz[i,it]
#         wfd2.vx[j] = bnd.vx[i,it]
#     end
#
#     # update pressure
#     A_mul_b!(tmp1, rfds.MpxBvx, wfd2.vx)
#     A_mul_b!(tmp2, rfds.MpzBvz, wfd2.vz)
#     addition!(tmp1, tmp2)
#     addition!(wfd2.p, wfd1.p, tmp1)
#
#     # correct for boundary part
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd2.p[j] = bnd.p[i,it]
#     end
#
#     return nothing
# end
#
# """
# save snapshot or full wave field or pressure field to hard drive, inject single source
# save_type="snapshot" : vz, vx, pz, px include PML part,
# save_type="wavefield": vz, vx, pz+px  without PML part,
# save_type="pressure" : p=pz+px without PML bounary part.
# """
# function multi_step_forward(path::String, src::Source, ofds::ObsorbFDStencil,
#          params::TdParams; save_type="pressure")
#
#     # initialize some variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp1 = zeros(params.data_format, params.Nz * params.Nx)
#     tmp2 = zeros(params.data_format, params.Nz * params.Nx)
#
#     # add source to the first snapshot
#     add_source!(spt1, src, 1)
#
#     # save snapshot
#     if save_type == "snapshot"
#        hdr = snapshot_header(params)
#        fid = write_RSheader(path, hdr)
#        append_one_snapshot(fid, spt1)
#
#     # save wavefield
#     elseif save_type == "wavefield"
#        hdr = wavefield_header(params)
#        fid = write_RSheader(path, hdr)
#        append_one_wavefield(fid, spt1, params)
#
#     # save pressure
#     elseif save_type == "pressure"
#        hdr = pressure_header(params)
#        fid = write_RSheader(path, hdr)
#        append_one_pressure(fid, spt1, params)
#     else
#        error("save_type only support snapshot, wavefield or pressure")
#     end
#
#     # loop over time-stepping
#     for it = 2 : params.nt
#
#         one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
#         add_source!(spt2, src, it)
#         copy_snapshot!(spt1, spt2)
#
#         if save_type == "snapshot"
#            append_one_snapshot(fid, spt1)
#         elseif save_type == "wavefield"
#            append_one_wavefield(fid, spt1, params)
#         elseif save_type == "pressure"
#            append_one_pressure(fid, spt1, params)
#         end
#     end
#
#     close(fid)
#     return nothing
# end
#
# """
#    compute recordings with a single source
# """
# function multi_step_forward!(rec::Recordings, src::Source,
#          ofds::ObsorbFDStencil, params::TdParams)
#
#     # initialize intermediate variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp1 = zeros(params.data_format, params.Nz * params.Nx)
#     tmp2 = zeros(params.data_format, params.Nz * params.Nx)
#
#     # add source to the first snapshot
#     add_source!(spt1, src, 1)
#
#     # obtain the first time sample of recordings
#     sample_spt2rec!(rec, spt1, 1)
#
#     # loop over time stepping
#     for it = 2 : params.nt
#
#         one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
#         add_source!(spt2, src, it)
#         sample_spt2rec!(rec, spt2, it)
#
#         # prepare for next time step
#         copy_snapshot!(spt1, spt2)
#     end
#
#     return nothing
# end
#
# """
#    compute recordings with multiple simultaneous source, the only difference is that the input
# is a vector of source
# """
# function multi_step_forward!(rec::Recordings, srcs::Vector{Source},
#          ofds::ObsorbFDStencil, params::TdParams)
#
#     # initialize intermediate variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp1 = zeros(params.data_format, params.Nz * params.Nx)
#     tmp2 = zeros(params.data_format, params.Nz * params.Nx)
#
#     # add multi-sources to the first snapshot
#     add_multi_sources!(spt1, srcs, 1)
#
#     # sample the snapshot to fill recordings
#     sample_spt2rec!(rec, spt1, 1)
#
#     # loop over time
#     for it = 2 : params.nt
#
#         one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
#         add_multi_sources!(spt2, srcs, it)
#         sample_spt2rec!(rec, spt2, it)
#
#         # prepare for the next time step
#         copy_snapshot!(spt1, spt2)
#     end
#
#     return nothing
# end
#
# """
#    compute wavefield boundary at each time step and the wavefield at last time step
# save them to hard drive.
# """
# function get_boundary_wavefield(path_bnd::String, path_wfd::String,
#          src::Source, ofds::ObsorbFDStencil, params::TdParams)
#
#     # one extra time step
#     nt = params.nt + 1
#
#     # initialize variables for time stepping
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp1 = zeros(params.data_format, params.Nz * params.Nx)
#     tmp2 = zeros(params.data_format, params.Nz * params.Nx)
#
#     # add source to the first snapshot
#     add_source!(spt1, src, 1)
#
#     # save the boundary of the first wavefield
#     hdr = boundary_header(params)
#     fid = write_RSheader(path_bnd, hdr)
#     append_one_boundary(fid, spt1, params)
#
#     # loop over time stepping
#     for it = 2 : nt
#
#         one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
#         add_source!(spt2, src, it)
#
#         # save boundary of wavefield
#         append_one_boundary(fid, spt2, params)
#
#         # prepare for next iteration
#         copy_snapshot!(spt1, spt2)
#     end
#
#     # finish saving all the boundary values, close the file
#     close(fid)
#
#     # save the wavefield at the final time step
#     wfd = sample_spt2wfd(spt1, params)
#     write_wavefield(path_wfd, wfd, params)
#
#     return nothing
# end
#
# """
#    Forward reconstruct acoustic pressure field using the saved boundary wavefield value and
# the source (used for concept proofing)
# """
# function pressure_reconstruct_forward(bnd::WavefieldBound, rfds::RigidFDStencil,
#          src::Source, params::TdParams)
#
#     # length of one-step pressure field
#     N = params.nz * params.nx
#
#     # decide the number of boundary layers
#     order = length(params.fd_coefficients)
#
#     # initialize intermediate variables
#     wfd1 = Wavefield(params)
#     wfd2 = Wavefield(params)
#     tmp1 = zeros(params.data_format, N)
#     tmp2 = zeros(params.data_format, N)
#
#     # allocate memory for saving pressure field
#     pre = zeros(params.data_format, N * params.nt)
#
#     # add source to the first wavefield
#     add_source!(wfd1, src, 1)
#
#     # save pressure field
#     copyto!(pre, 1, wfd1.p, 1, N)
#
#     # the start index for saving the next pressure field
#     idx_o = N+1
#
#     # loop over time stepping
#     for it = 2 : params.nt
#
#         # forward time steping and correcting boundaries
#         one_step_forward!(wfd2, wfd1, rfds, it, bnd, tmp1, tmp2, params)
#
#         # add source to wavefield
#         if order+1 <= src.isz <= params.nz-order && order+1 <= src.isx <= params.nx-order
#            add_source!(wfd2, src, it)
#         end
#
#         # save the pressure field
#         copyto!(pre, idx_o, wfd2.p, 1, N)
#         idx_o = idx_o + N
#
#         # prepare for next step
#         copy_wavefield!(wfd1, wfd2)
#     end
#
#     return reshape(pre, params.nz, params.nx, params.nt)
# end

# completely based on for loop
# function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, params::TdParams{Ti,Tv},
#          fc::Vector{Tv}, tmp::Vector{Tv}, tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
#          tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # # local variables
#     # fc = params.fc
#
#     # number of byte for each element
#     esize = sizeof(params.data_format)
#
#     # save the pressure field
#     tmp .= spt1.pz .+ spt1.px
#     p    = pointer(tmp)
#
#     # update vz=================================================================
#     ilower = 1
#     iupper = params.Nz
#     for ix = 1 : params.Nx
#
#         # copy one column
#         copyto!(tmp_z1, 1, tmp, ilower, params.Nz)
#         fill!(tmp_z2, 0.0)
#
#         # dpdz
#         for iz = 1 : params.order-1
#             for idx = 1 : iz
#                 tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx+1]
#             end
#             for idx = 1 : params.order
#                 tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx]
#             end
#         end
#
#         for iz = params.order : params.Nz-params.order
#             for idx = 1 : params.order
#                 tmp_z2[iz] += fc[idx] * (tmp_z1[iz+idx] - tmp_z1[iz-idx+1])
#             end
#         end
#
#         for iz = params.Nz-params.order+1 : params.Nz
#             for idx = 1 : params.order
#                 tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx+1]
#             end
#             for idx = 1 : params.Nz - iz
#                 tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx]
#             end
#         end
#
#         # update one column
#         idx    = 0
#         for iz = ilower : iupper
#             idx = idx + 1
#             spt2.vz[iz] = params.MvzBvz[idx]*spt1.vz[iz] + params.MvzBp[iz]*tmp_z2[idx]
#         end
#
#         # prepare for next column
#         ilower = ilower + params.Nz
#         iupper = iupper + params.Nz
#     end
#
#     # update vx=================================================================
#     shift = 0
#     for iz = 1 : params.Nz
#
#         # copy one row
#         BLAS.blascopy!(params.Nx, p+shift, params.Nz, tmp_x1, 1)
#         fill!(tmp_x2, 0.0)
#
#         # dpdx
#         for ix = 1 : params.order-1
#             for idx = 1 : ix
#                 tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx+1]
#             end
#             for idx = 1 : params.order
#                 tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx]
#             end
#         end
#
#         for ix = params.order : params.Nx-params.order
#             for idx = 1 : params.order
#                 tmp_x2[ix] += fc[idx] * (tmp_x1[ix+idx] - tmp_x1[ix-idx+1])
#             end
#         end
#
#         for ix = params.Nx-params.order+1 : params.Nx
#             for idx = 1 : params.order
#                 tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx+1]
#             end
#             for idx = 1 : params.Nx - ix
#                 tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx]
#             end
#         end
#
#         idx = iz
#         for ix = 1 : params.Nx
#             spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*tmp_x2[ix]
#             idx= idx + params.Nz
#         end
#
#         # prepare for next row
#         shift = shift + esize
#     end
#
#     # update pz ================================================================
#     ilower = 1
#     iupper = params.Nz
#
#     # loop over columns
#     for ix = 1 : params.Nx
#
#         # copy one column
#         copyto!(tmp_z1, 1, spt2.vz, ilower, params.Nz)
#         fill!(tmp_z2, 0.0)
#
#         # dvdz
#         for iz = 1 : params.order
#             for idx = 1 : iz-1
#                 tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx]
#             end
#             for idx = 1 : params.order
#                 tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx-1]
#             end
#         end
#
#         for iz = params.order+1 : params.Nz-params.order+1
#             for idx = 1 : params.order
#                 tmp_z2[iz] += fc[idx] * (tmp_z1[iz+idx-1]-tmp_z1[iz-idx])
#             end
#         end
#
#         for iz = param.Nz-params.order+2 : params.Nz
#             for idx = 1 : params.order
#                 tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx]
#             end
#             for idx = 1 : params.Nz-iz+1
#                 tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx-1]
#             end
#         end
#
#         idx    = 0
#         for iz = ilower : iupper
#             idx = idx + 1
#             spt2.pz[iz] = params.MpzBpz[idx]*spt1.pz[iz] + params.MpzBvz[iz]*tmp_z2[idx]
#         end
#
#         # prepare for next column
#         ilower = ilower + params.Nz
#         iupper = iupper + params.Nz
#     end
#
#     # update px=================================================================
#     shift = 0
#     p = pointer(spt2.vx)
#     for iz = 1 : params.Nz
#
#         # copy one row
#         BLAS.blascopy!(params.Nx, p+shift, params.Nz, tmp_x1, 1)
#         fill!(tmp_x2, 0.0)
#
#         # dvdx
#         for ix = 1 : params.order
#             for idx = 1 : ix
#                 tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx]
#             end
#             for idx = 1 : params.order
#                 tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx-1]
#             end
#         end
#
#         for ix = params.order+1 : params.Nx-params.order+1
#             for idx = 1 : params.order
#                 tmp_x2[ix] += fc[idx] * (tmp_x1[ix+idx-1] - tmp_x1[ix-idx])
#             end
#         end
#
#         for ix = params.Nx-params.order+2 : params.Nx
#             for idx = 1 : params.order
#                 tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx]
#             end
#             for idx = 1 : params.Nx - ix + 1
#                 tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx-1]
#             end
#         end
#
#         idx = iz
#         for ix = 1 : params.Nx
#             spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*tmp_x2[ix]
#             idx= idx + params.Nz
#         end
#
#         # prepare for next row
#         shift = shift + esize
#     end
#
# end
