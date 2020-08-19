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
function multi_step_adjoint(rec::Recordings, params::TdParams;
                             path_spt="NULL", path_wfd="NULL" , path_pre="NULL", interval=1)

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp    = zeros(params.data_format, params.Nz * params.Nx)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # inject the recordings to the last snapshot
    inject_rec2spt!(spt2, rec, params.nt)

    # the time range of source function
    if path_spt != "NULL"
       hdr_spt = snapshot_header(params, interval)
       fid_spt = write_RSheader(path_spt, hdr_spt)
       append_one_snapshot(fid_spt, spt2)
    end

    # save wavefield
    if path_wfd != "NULL"
       hdr_wfd = wavefield_header(params, interval)
       fid_wfd = write_RSheader(path_wfd, hdr_wfd)
       append_one_wavefield(fid_wfd, spt2, params)
    end

    # save pressure
    if path_pre != "NULL"
       hdr_pre = pressure_header(params, interval)
       fid_pre = write_RSheader(path_pre, hdr_pre)
       append_one_pressure(fid_pre, spt2, params)
    end

    # back propagation
    for it = params.nt-1 : -1 : 1

        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, rec, it)
        copy_snapshot!(spt2, spt1)

        # save wavefield at current step to disk
        if mod(it-1, interval) == 0
           path_spt != "NULL" && append_one_snapshot(fid_spt, spt2)
           path_wfd != "NULL" && append_one_wavefield(fid_wfd, spt2, params)
           path_pre != "NULL" && append_one_pressure(fid_pre, spt2, params)
        end
    end

    # close file and flush the buffer
    path_spt != "NULL" && close(fid_spt)
    path_wfd != "NULL" && close(fid_wfd)
    path_pre != "NULL" && close(fid_pre)

    # reverse the time order of the adjoint wavefield
    if path_spt != "NULL"
       tmp_dir = splitdir(path_spt)
       path_tmp= joinpath(tmp_dir[1], "spt_reverse_file.rsb")
       reverse_order(path_tmp, path_spt; save_flag="snapshot")
       mv(path_tmp, path_spt, force=true)
    end

    if path_wfd != "NULL"
       tmp_dir = splitdir(path_wfd)
       path_tmp= joinpath(tmp_dir[1], "wfd_reverse_file.rsb")
       reverse_order(path_tmp, path_wfd; save_flag="wavefield")
       mv(path_tmp, path_wfd, force=true)
    end

    if path_pre != "NULL"
       tmp_dir = splitdir(path_pre)
       path_tmp= joinpath(tmp_dir[1], "pre_reverse_file.rsb")
       reverse_order(path_tmp, path_pre; save_flag="pressure")
       mv(path_tmp, path_pre, force=true)
    end

    return nothing
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
function pressure_reconstruct_backward(path_bnd::Tp, path_lwfd::Tp,
         src::Ts, params::TdParams) where {Tp<:String, Ts<:Union{Source,Vector{Source}}}

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
    wfd2 = read_one_wavefield(path_lwfd, 1)
    tmp_z1 = zeros(params.data_format, params.nz)
    tmp_z2 = zeros(params.data_format, params.nz)
    tmp_x1 = zeros(params.data_format, params.nx)
    tmp_x2 = zeros(params.data_format, params.nx)

    # initialize boundary value
    fid_bnd = open(path_bnd, "r")
    bnd     = WavefieldBound(params)

    # pre-allocate memory for saving the reconstructed pressure field
    pre = zeros(params.data_format, N * params.nt)
    idx_o = N*params.nt - N + 1 # the pressure field at last time step
    copyto!(pre, idx_o, wfd2.p, 1, N)

    # prepare for backward reconstruction
    subtract_source!(wfd2, srcs, params.nt)

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
        subtract_source!(wfd2, srcs, it)
    end

    close(fid_bnd)
    return reshape(pre, params.nz, params.nx, params.nt)
end
