"""
   The adjoint operator of one_step_forward (in-place).
"""
function one_step_adjoint!(spt2::Snapshot, spt1::Snapshot, params::TdParams,
         tmp::Vector{Tv}, tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
         tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Tv<:AbstractFloat}

    # update vz column-by-column
    for ix = 1 : params.Nx
        ilower = (ix-1) * params.Nz + 1
        iupper = ilower + params.Nz - 1

        idx = 0
        for iz = ilower : iupper
            idx = idx + 1
            tmp_z1[idx] = params.MpzBvz[iz] * spt1.pz[iz]
        end

        At_mul_b!(tmp_z2, params.dvdz, tmp_z1)
        # mul!(tmp_z2, transpose(params.dvdz), tmp_z1)

        idx = 0
        for iz = ilower : iupper
            idx = idx + 1
            spt2.vz[iz] = spt1.vz[iz] + tmp_z2[idx]
            spt2.pz[iz] = spt1.pz[iz] * params.MpzBpz[idx]
        end
    end

    # update vx row-by-row
    for iz = 1 : params.Nz

        idx = iz
        for ix = 1 : params.Nx
            tmp_x1[ix] = params.MpxBvx[idx] * spt1.px[idx]
            idx= idx + params.Nz
        end

        At_mul_b!(tmp_x2, params.dvdx, tmp_x1)
        # mul!(tmp_x2, transpose(params.dvdx), tmp_x1)

        idx = iz
        for ix = 1 : params.Nx
            spt2.vx[idx] = spt1.vx[idx] + tmp_x2[ix]
            spt2.px[idx] = spt1.px[idx] * params.MpxBpx[ix]
            idx= idx + params.Nz
        end
    end

    # update pz column-by-column
    for ix = 1 : params.Nx
        ilower = (ix-1) * params.Nz + 1
        iupper = ilower + params.Nz - 1

        idx = 0
        for iz = ilower : iupper
            idx = idx + 1
            tmp_z1[idx] = params.MvzBp[iz] * spt2.vz[iz]
        end

        At_mul_b!(tmp_z2, params.dpdz, tmp_z1)
        # mul!(tmp_z2, transpose(params.dpdz), tmp_z1)

        idx = 0
        for iz = ilower : iupper
            idx = idx + 1
            tmp[iz] = tmp_z2[idx]
        end
    end

    # update px row-by-row
    for iz = 1 : params.Nz

        idx = iz
        for ix = 1 : params.Nx
            tmp_x1[ix] = params.MvxBp[idx] * spt2.vx[idx]
            idx= idx + params.Nz
        end

        At_mul_b!(tmp_x2, params.dpdx, tmp_x1)
        # mul!(tmp_x2, transpose(params.dpdx), tmp_x1)

        idx = iz
        for ix = 1 : params.Nx
            tmp[idx] = tmp[idx] + tmp_x2[ix]
            idx= idx + params.Nz
        end
    end

    # finalize the result
    for ix = 1 : params.Nx
        row_idx = (ix-1) * params.Nz

        for iz = 1 : params.Nz
            idx = row_idx + iz

            spt2.vz[idx] = params.MvzBvz[iz] * spt2.vz[idx]
            spt2.vx[idx] = params.MvxBvx[ix] * spt2.vx[idx]
            spt2.pz[idx] = tmp[idx] + spt2.pz[idx]
            spt2.px[idx] = tmp[idx] + spt2.px[idx]
        end
    end

    return nothing
end

"""
   save the adjoint snapshot, wavefield or pressure field
"""
function multi_step_adjoint!(path::String, rec::Recordings, params::TdParams;
                             save_flag="pressure")

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp  = zeros(params.data_format, params.Nz * params.Nx)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # inject the recordings to the last snapshot
    inject_rec2spt!(spt2, rec, params.nt)

    # the time range of source function
    if save_flag == "snapshot"
       hdr = snapshot_header(params)
       fid = write_RSheader(path, hdr)
       append_one_snapshot(fid, spt2)

    # save wavefield
    elseif save_flag == "wavefield"
       hdr = wavefield_header(params)
       fid = write_RSheader(path, hdr)
       append_one_wavefield(fid, spt2, params)

    # save pressure
    elseif save_flag == "pressure"
       hdr = pressure_header(params)
       fid = write_RSheader(path, hdr)
       append_one_pressure(fid, spt2, params)
    else
       error("save_type only support snapshot, wavefield or pressure")
    end

    # back propagation
    for it = params.nt-1 : -1 : 1

        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, rec, it)
        copy_snapshot!(spt2, spt1)

        if save_flag == "snapshot"
           append_one_snapshot(fid, spt2)

        elseif save_flag == "wavefield"
           append_one_wavefield(fid, spt2, params)

        elseif save_flag == "pressure"
           append_one_pressure(fid, spt2, params)
        end
    end

    close(fid)

    # reverse the time order of the adjoint wavefield
    tmp_dir = splitdir(path)
    path_tmp= joinpath(tmp_dir[1], "tmp_reverse_file.rsb")
    reverse_order(path_tmp, path; save_flag=save_flag)
    mv(path_tmp, path, force=true)

    return nothing
end

"""
   the adjoint operator of one point source simulation, only used for testing the
dot product test
"""
function multi_step_adjoint!(rec::Recordings, src::Source, params::TdParams)

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp  = zeros(params.data_format, params.Nz * params.Nx)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # inject the recordings to the last snapshot
    inject_rec2spt!(spt2, rec, params.nt)

    # the time range of source function
    p   = zeros(params.data_format, src.it_max - src.it_min + 1)
    idx = src.src2spt

    if src.it_min <= params.nt <= src.it_max
       i = params.nt - src.it_min + 1
       p[i] = spt2.px[idx] * src.dt
    end

    # back propagation
    for it = params.nt-1 : -1 : 1

        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, rec, it)

        if src.it_min <= it <= src.it_max
           i = it - src.it_min + 1
           p[i] = spt1.px[idx] * src.dt
        end

        copy_snapshot!(spt2, spt1)
    end

    return p
end

"""
   one step backward reconstruction of wave field from the boundary values, the boundary
wavefield at only one time step is provided
"""
function one_step_backward!(wfd2::Wavefield, wfd1::Wavefield,
         bnd::WavefieldBound, params::TdParams,
         tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
         tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    # update pressure column-by-column
    for ix = 1 : params.nx
        ilower = (ix-1) * params.nz + 1
        iupper = ilower + params.nz - 1

        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = wfd1.vz[iz]
        end

        A_mul_b!(tmp_z2, params.rvdz, tmp_z1)        # dvdz

        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            wfd2.p[iz] = wfd1.p[iz] - params.RpBv[iz]*tmp_z2[idx]
        end
    end

    # update pressure row-by-row
    for iz = 1 : params.nz

        idx = iz
        for ix = 1 : params.nx
            tmp_x1[ix] = wfd1.vx[idx]
            idx= idx + params.nz
        end

        A_mul_b!(tmp_x2, params.rvdx, tmp_x1)         # dvdx

        idx = iz
        for ix = 1 : params.nx
            wfd2.p[idx] = wfd2.p[idx] - params.RpBv[idx]*tmp_x2[ix]
            idx= idx + params.nz
        end
    end

    # correct for the boundary value
    for i = 1 : length(params.wfd2bnd)
        j = params.wfd2bnd[i]
        wfd2.p[j] = bnd.p[i]
    end

    # update vz column-by-column
    for ix = 1 : params.nx
        ilower = (ix-1) * params.nz + 1
        iupper = ilower + params.nz - 1

        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = wfd2.p[iz]
        end

        A_mul_b!(tmp_z2, params.rpdz, tmp_z1)        # dpdz

        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            wfd2.vz[iz] = wfd1.vz[iz] - params.RvzBp[iz]*tmp_z2[idx]
        end
    end

    # update vx row-by-row
    for iz = 1 : params.nz

        idx = iz
        for ix = 1 : params.nx
            tmp_x1[ix] = wfd2.p[idx]
            idx= idx + params.nz
        end

        A_mul_b!(tmp_x2, params.rpdx, tmp_x1)         # dpdx

        idx = iz
        for ix = 1 : params.nx
            wfd2.vx[idx] = wfd1.vx[idx] - params.RvxBp[idx]*tmp_x2[ix]
            idx= idx + params.nz
        end
    end

    # correct for boundary part
    for i = 1 : length(params.wfd2bnd)
        j = params.wfd2bnd[i]
        wfd2.vz[j] = bnd.vz[i]
        wfd2.vx[j] = bnd.vx[i]
    end

    return nothing
end

"""
   Reconstruct pressure field backward using the boundary wavefield value and
the last wavefield
"""
function pressure_reconstruct_backward(path_bnd::Ts, path_wfd::Ts,
         src::Source, params::TdParams) where {Ts <: String}

    # length of one-step pressure field
    N = params.nz * params.nx

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = read_one_wavefield(path_wfd, 1)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd = WavefieldBound(params)

    # pre-allocate memory for saving the reconstructed pressure field
    pre = zeros(params.data_format, N * params.nt)
    idx_o = N*params.nt - N + 1 # the pressure field at last time step
    copyto!(pre, idx_o, wfd2.p, 1, N)

    # prepare for backward reconstruction
    subtract_source!(wfd2, src, params.nt)

    # backward iterations
    for it = params.nt-1 : -1 : 1

        # read the boundary value
        read_one_boundary!(bnd, fid_bnd, it, params)

        # one step backward reconstruction
        one_step_backward!(wfd1, wfd2, bnd, params,
                           tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # save the pressure field
        idx_o = idx_o - N
        copyto!(pre, idx_o, wfd1.p, 1, N)

        # prepare for next step
        copy_wavefield!(wfd2, wfd1)
        subtract_source!(wfd2, src, it)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt)
end

"""
   Reconstruct pressure field backward using the boundary wavefield value and
the last wavefield when the wavefield is generated with multiple sources
"""
function pressure_reconstruct_backward(path_bnd::Ts, path_wfd::Ts,
         srcs::Vector{Source}, params::TdParams) where {Ts <: String}

    # length of one-step pressure field
    N = params.nz * params.nx

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = read_one_wavefield(path_wfd, 1)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd = WavefieldBound(params)

    # pre-allocate memory for saving the reconstructed pressure field
    pre = zeros(params.data_format, N * params.nt)
    idx_o = N*params.nt - N + 1                     # the pressure field at last time step
    copyto!(pre, idx_o, wfd2.p, 1, N)

    # prepare for the backward reconstruction
    subtract_multi_sources!(wfd2, srcs, params.nt)

    # backward iterations
    for it = params.nt-1 : -1 : 1

        # read the boundary value
        read_one_boundary!(bnd, fid_bnd, it, params)

        # one step back
        one_step_backward!(wfd1, wfd2, bnd, params,
                           tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # save the pressure field
        idx_o = idx_o - N
        copyto!(pre, idx_o, wfd1.p, 1, N)

        # prepare for next step
        copy_wavefield!(wfd2, wfd1)
        subtract_multi_sources!(wfd2, srcs, it)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt)
end

# """
#    one step backward reconstruction of wave field from the boundary values
# """
# function one_step_backward!(wfd2::Wavefield, wfd1::Wavefield, it::Ti,
#          bnd::WavefieldBound, params::TdParams,
#          tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
#          tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # update pressure column-by-column
#     for ix = 1 : params.nx
#         ilower = (ix-1) * params.nz + 1
#         iupper = ilower + params.nz - 1
#
#         idx    = 0
#         for iz = ilower : iupper
#             idx= idx + 1
#             tmp_z1[idx] = wfd1.vz[iz]
#         end
#
#         A_mul_b!(tmp_z2, params.rvdz, tmp_z1)        # dvdz
#
#         idx    = 0
#         for iz = ilower : iupper
#             idx = idx + 1
#             wfd2.p[iz] = wfd1.p[iz] - params.RpzBvz[iz]*tmp_z2[idx]
#         end
#     end
#
#     # update pressure row-by-row
#     for iz = 1 : params.nz
#
#         idx = iz
#         for ix = 1 : params.nx
#             tmp_x1[ix] = wfd1.vx[idx]
#             idx= idx + params.nz
#         end
#
#         A_mul_b!(tmp_x2, params.rvdx, tmp_x1)         # dvdx
#
#         idx = iz
#         for ix = 1 : params.nx
#             wfd2.p[idx] = wfd2.p[idx] - params.RpxBvx[idx]*tmp_x2[ix]
#             idx= idx + params.nz
#         end
#     end
#
#     # correct for the boundary value
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd2.p[j] = bnd.p[i,it]
#     end
#
#     # update vz column-by-column
#     for ix = 1 : params.nx
#         ilower = (ix-1) * params.nz + 1
#         iupper = ilower + params.nz - 1
#
#         idx    = 0
#         for iz = ilower : iupper
#             idx= idx + 1
#             tmp_z1[idx] = wfd2.p[iz]
#         end
#
#         A_mul_b!(tmp_z2, params.rpdz, tmp_z1)        # dpdz
#
#         idx    = 0
#         for iz = ilower : iupper
#             idx = idx + 1
#             wfd2.vz[iz] = wfd1.vz[iz] - params.RvzBp[iz]*tmp_z2[idx]
#         end
#     end
#
#     # update vx row-by-row
#     for iz = 1 : params.nz
#
#         idx = iz
#         for ix = 1 : params.nx
#             tmp_x1[ix] = wfd2.p[idx]
#             idx= idx + params.nz
#         end
#
#         A_mul_b!(tmp_x2, params.rpdx, tmp_x1)         # dpdx
#
#         idx = iz
#         for ix = 1 : params.nx
#             wfd2.vx[idx] = wfd1.vx[idx] - params.RvxBp[idx]*tmp_x2[ix]
#             idx= idx + params.nz
#         end
#     end
#
#     # correct for boundary part
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd2.vz[j] = bnd.vz[i,it]
#         wfd2.vx[j] = bnd.vx[i,it]
#     end
#
#     return nothing
# end

# """
#    Reconstruct pressure field backward using the boundary wavefield value and
# the last wavefield
# """
# function pressure_reconstruct_backward(bnd::WavefieldBound, wfd::Wavefield,
#          src::Source, params::TdParams)
#
#     # length of one-step pressure field
#     N = params.nz * params.nx
#
#     # initialize intermediate variables
#     wfd1 = Wavefield(params)
#     wfd2 = Wavefield(params)
#     tmp_z1 = zeros(params.data_format, params.nz)
#     tmp_z2 = zeros(params.data_format, params.nz)
#     tmp_x1 = zeros(params.data_format, params.nx)
#     tmp_x2 = zeros(params.data_format, params.nx)
#
#     # the last snapshot
#     copy_wavefield!(wfd2, wfd)
#
#     # pre-allocate memory for saving the reconstructed pressure field
#     pre = zeros(params.data_format, N * params.nt)
#     idx_o = N*params.nt - N + 1 # the pressure field at last time step
#     copyto!(pre, idx_o, wfd2.p, 1, N)
#
#     # prepare for backward reconstruction
#     subtract_source!(wfd2, src, params.nt)
#
#     # backward iterations
#     for it = params.nt-1 : -1 : 1
#
#         # one step back
#         one_step_backward!(wfd1, wfd2, it, bnd, params,
#                            tmp_z1, tmp_z2, tmp_x1, tmp_x2)
#
#         # save the pressure field
#         idx_o = idx_o - N
#         copyto!(pre, idx_o, wfd1.p, 1, N)
#
#         # prepare for next step
#         copy_wavefield!(wfd2, wfd1)
#         subtract_source!(wfd2, src, it)
#     end
#
#     return reshape(pre, params.nz, params.nx, params.nt)
# end

# """
#    The adjoint operator of one_step_forward.
# """
# function one_step_adjoint!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil)
#
#     spt2.vz = spt1.vz + (fidMtx.MpzBvz)' * spt1.pz
#     spt2.vx = spt1.vx + (fidMtx.MpxBvx)' * spt1.px
#     spt2.pz =            fidMtx.MpzBpz  .* spt1.pz
#     spt2.px =            fidMtx.MpxBpx  .* spt1.px
#
#     spt2.pz = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.pz
#     spt2.px = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.px
#     spt2.vz =  fidMtx.MvzBvz  .* spt2.vz
#     spt2.vx =  fidMtx.MvxBvx  .* spt2.vx
#
#     return nothing
# end
#
# """
#    The adjoint operator of one_step_forward (in-place).
# """
# function one_step_adjoint!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil,
#          tmp1::Vector{Tv}, tmp2::Vector{Tv}) where {Tv<:AbstractFloat}
#
#     At_mul_b!(tmp1, ofds.MpzBvz, spt1.pz)
#     addition!(spt2.vz, spt1.vz, tmp1)
#
#     At_mul_b!(tmp1, ofds.MpxBvx, spt1.px)
#     addition!(spt2.vx, spt1.vx, tmp1)
#
#     multiplication!(spt2.pz, ofds.MpzBpz, spt1.pz)
#     multiplication!(spt2.px, ofds.MpxBpx, spt1.px)
#
#     At_mul_b!(tmp1, ofds.MvzBp, spt2.vz)
#     At_mul_b!(tmp2, ofds.MvxBp, spt2.vx)
#     addition!(tmp1, tmp2)
#
#     addition!(spt2.pz, tmp1)
#     addition!(spt2.px, tmp1)
#
#     multiplication!(spt2.vz, ofds.MvzBvz)
#     multiplication!(spt2.vx, ofds.MvxBvx)
#
#     return nothing
# end
#
# """
#    one step backward reconstruction of wave field from the boundary values
# """
# function one_step_backward!(wfd1::Wavefield, wfd2::Wavefield, rfds::RigidFDStencil, it::Ti,
#          bnd::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, params::TdParams) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # wfd1.p  = wfd2.p  - rfds.MpxBvx * wfd2.vx - rfds.MpzBvz * wfd2.vz
#     # wfd1.p[bound.index]  = bound.p[:, it]
#     A_mul_b!(tmp1, rfds.MpxBvx, wfd2.vx)
#     A_mul_b!(tmp2, rfds.MpzBvz, wfd2.vz)
#     addition!(tmp1, tmp2)
#     minus!(wfd1.p, wfd2.p, tmp1)
#
#     # correct for boundary part of pressure
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd1.p[j] = bnd.p[i,it]
#     end
#
#     # wfd1.vx = wfd2.vx - rfds.MvxBp  * wfd1.p
#     # wfd1.vz = wfd2.vz - rfds.MvzBp  * wfd1.p
#     # wfd1.vx[bound.index] = bound.vx[:, it]
#     # wfd1.vz[bound.index] = bound.vz[:, it]
#     A_mul_b!(tmp1, rfds.MvxBp, wfd1.p)
#     A_mul_b!(tmp2, rfds.MvzBp, wfd1.p)
#     minus!(wfd1.vx, wfd2.vx, tmp1)
#     minus!(wfd1.vz, wfd2.vz, tmp2)
#
#     # correct for boundary part of particle velocity
#     for i = 1 : length(params.wfd2bnd)
#         j = params.wfd2bnd[i]
#         wfd1.vz[j] = bnd.vz[i,it]
#         wfd1.vx[j] = bnd.vx[i,it]
#     end
#
#     return nothing
# end
#
# """
#    the adjoint operator of one point source simulation.
# """
# function multi_step_adjoint(rec::Recordings, ofds::ObsorbFDStencil,
#          src::Source, params::TdParams)
#
#     # initialize intermediate variables
#     spt1 = Snapshot(params)
#     spt2 = Snapshot(params)
#     tmp1 = zeros(params.data_format, params.Nz * params.Nx)
#     tmp2 = zeros(params.data_format, params.Nz * params.Nx)
#
#     # inject the recordings to the last snapshot
#     inject_rec2spt!(spt2, rec, params.nt)
#
#     # the time range of source function
#     p   = zeros(params.data_format, src.it_max - src.it_min + 1)
#     idx = src.src2spt
#
#     if src.it_min <= params.nt <= src.it_max
#        i = params.nt - src.it_min + 1
#        p[i] = spt2.pz[idx] * params.dt
#     end
#
#     # back propagation
#     for it = params.nt-1 : -1 : 1
#
#         one_step_adjoint!(spt1, spt2, ofds, tmp1, tmp2)
#         inject_rec2spt!(spt1, rec, it)
#
#         if src.it_min <= it <= src.it_max
#            i = it - src.it_min + 1
#            p[i] = spt1.pz[idx] * params.dt
#         end
#
#         copy_snapshot!(spt2, spt1)
#     end
#
#     return p
# end
#
# """
#    Reconstruct pressure field backward using the boundary wavefield value and
# the last wavefield
# """
# function pressure_reconstruct_backward(bnd::WavefieldBound, wfd::Wavefield,
#          rfds::RigidFDStencil, src::Source, params::TdParams)
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
#     copy_wavefield!(wfd2, wfd)
#
#     # pre-allocate memory for saving the reconstructed pressure field
#     pre = zeros(params.data_format, N * params.nt)
#     idx_o = N*params.nt - N + 1 # the pressure field at last time step
#
#     for it = params.nt : -1 : 1
#
#         # one step back
#         one_step_backward!(wfd1, wfd2, rfds, it, bnd, tmp1, tmp2, params)
#
#         # save the pressure field
#         copyto!(pre, idx_o, wfd1.p, 1, N)
#         idx_o = idx_o - N
#
#         # prepare for next step
#         copy_wavefield!(wfd2, wfd1)
#         subtract_source!(wfd2, src, it)
#     end
#
#     return reshape(pre, params.nz, params.nx, params.nt)
# end
