function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, params::ModelParams{Ti,Tv},
         tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
         tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    # update vz column-by-column
    for ix = 1 : params.Nx
        ilower = (ix-1) * params.Nz + 1
        iupper = ilower + params.Nz - 1
        idx    = 0
        for iz = ilower : iupper
            idx= idx + 1
            tmp_z1[idx] = spt1.pz[iz] + spt1.px[iz]  # p1 = pz1 + px1
        end
        A_mul_b!(tmp_z2, params.dpdz, tmp_z1)        # dpdz
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
        A_mul_b!(tmp_x2, params.dpdx, tmp_x1)         # dpdx
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
        idx = iz
        for ix = 1 : params.Nx
            spt2.px[idx] = params.MpxBpx[ix]*spt1.px[idx] + params.MpxBvx[idx]*tmp_x2[ix]
            idx= idx + params.Nz
        end
    end

    return nothing
end

"""
save snapshot or full wave field or pressure field to hard drive, inject single source
save_type="snapshot" : vz, vx, pz, px include PML part,
save_type="wavefield": vz, vx, pz+px  without PML part,
save_type="pressure" : p=pz+px without PML bounary part.
"""
function multi_step_forward!(path::String, src::Source, params::ModelParams;
                            save_flag="pressure")

    # initialize some variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # save snapshot
    if save_flag == "snapshot"
       hdr = snapshot_header(params)
       fid = write_RSheader(path, hdr)
       append_one_snapshot(fid, spt1)

    # save wavefield
    elseif save_flag == "wavefield"
       hdr = wavefield_header(params)
       fid = write_RSheader(path, hdr)
       append_one_wavefield(fid, spt1, params)

    # save pressure
    elseif save_flag == "pressure"
       hdr = pressure_header(params)
       fid = write_RSheader(path, hdr)
       append_one_pressure(fid, spt1, params)
    else
       error("save_type only support snapshot, wavefield or pressure")
    end

    # loop over time-stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        add_source!(spt2, src, it)
        copy_snapshot!(spt1, spt2)

        if save_flag == "snapshot"
           append_one_snapshot(fid, spt1)
        elseif save_flag == "wavefield"
           append_one_wavefield(fid, spt1, params)
        elseif save_flag == "pressure"
           append_one_pressure(fid, spt1, params)
        end
    end

    close(fid)
    return nothing
end

"""
save snapshot or full wave field or pressure field to hard drive, inject multiple simultaneouse source
save_type="snapshot" : vz, vx, pz, px include PML part,
save_type="wavefield": vz, vx, pz+px  without PML part,
save_type="pressure" : p=pz+px without PML bounary part.
"""
function multi_step_forward!(path::String, srcs::Vector{Source}, params::ModelParams;
                            save_flag="pressure")

    # initialize some variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_multi_sources!(spt1, srcs, 1)

    # save snapshot
    if save_flag == "snapshot"
       hdr = snapshot_header(params)
       fid = write_RSheader(path, hdr)
       append_one_snapshot(fid, spt1)

    # save wavefield
    elseif save_flag == "wavefield"
       hdr = wavefield_header(params)
       fid = write_RSheader(path, hdr)
       append_one_wavefield(fid, spt1, params)

    # save pressure
    elseif save_flag == "pressure"
       hdr = pressure_header(params)
       fid = write_RSheader(path, hdr)
       append_one_pressure(fid, spt1, params)
    else
       error("save_type only support snapshot, wavefield or pressure")
    end

    # loop over time-stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        add_multi_sources!(spt2, srcs, it)
        copy_snapshot!(spt1, spt2)

        if save_flag == "snapshot"
           append_one_snapshot(fid, spt1)
        elseif save_flag == "wavefield"
           append_one_wavefield(fid, spt1, params)
        elseif save_flag == "pressure"
           append_one_pressure(fid, spt1, params)
        end
    end

    close(fid)
    return nothing
end

"""
   compute the recordings via forward modelling with single source, the boundary and the wavefield
at the last time step are saved if "path_bnd" and "path_wfd" are given.
"""
function multi_step_forward!(rec::Recordings, src::Source, params::ModelParams;
         path_bnd="NULL", path_wfd="NULL")

    # initialize variables for time stepping
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # obtain the first time sample of recordings
    sample_spt2rec!(rec, spt1, 1)

    # save the boundary of the first wavefield
    if path_bnd != "NULL"
       hdr = boundary_header(params)
       fid = write_RSheader(path_bnd, hdr)
       append_one_boundary(fid, spt1, params)
    end

    # loop over time stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        add_source!(spt2, src, it)
        sample_spt2rec!(rec, spt1, it)

        # save boundary of wavefield
        if path_bnd != "NULL"
           append_one_boundary(fid, spt2, params)
        end

        # prepare for next iteration
        copy_snapshot!(spt1, spt2)
    end

    # finish saving all the boundary values, close the file
    if path_bnd != "NULL"
       close(fid)
    end

    # save the wavefield at the final time step
    if path_wfd != "NULL"
       wfd = sample_spt2wfd(spt1, params)
       write_wavefield(path_wfd, wfd, params)
    end

    return nothing
end

"""
   compute the recordings via forward modelling with simultaneouse source, the boundary and the wavefield
at the last time step are saved if "path_bnd" and "path_wfd" are given.
"""
function multi_step_forward!(rec::Recordings, srcs::Vector{Source}, params::ModelParams;
         path_bnd="NULL", path_wfd="NULL")

    # initialize variables for time stepping
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_multi_sources!(spt1, srcs, 1)

    # obtain the first time sample of recordings
    sample_spt2rec!(rec, spt1, 1)


    # save the boundary of the first wavefield
    if path_bnd != "NULL"
       hdr = boundary_header(params)
       fid = write_RSheader(path_bnd, hdr)
       append_one_boundary(fid, spt1, params)
    end

    # loop over time stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        add_multi_sources!(spt2, srcs, it)
        sample_spt2rec!(rec, spt1, it)

        # save boundary of wavefield
        if path_bnd != "NULL"
           append_one_boundary(fid, spt2, params)
        end

        # prepare for next iteration
        copy_snapshot!(spt1, spt2)
    end

    # finish saving all the boundary values, close the file
    if path_bnd != "NULL"
       close(fid)
    end

    # save the wavefield at the final time step
    if path_wfd != "NULL"
       wfd = sample_spt2wfd(spt1, params)
       write_wavefield(path_wfd, wfd, params)
    end

    return nothing
end

"""
   one step-forward reconstruction of wavefield via recording the boundary value
"""
function one_step_forward!(wfd2::Wavefield, wfd1::Wavefield,
         bnd::WavefieldBound, params::ModelParams,
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
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
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
            wfd2.p[iz] = wfd1.p[iz] + params.RpzBvz[iz]*tmp_z2[idx]
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
            wfd2.p[idx] = wfd2.p[idx] + params.RpxBvx[idx]*tmp_x2[ix]
            idx= idx + params.nz
        end
    end

    # correct for boundary part
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd2.p[j] = bnd.p[i]
    end

    return nothing
end

"""
   Forward reconstruct acoustic pressure field using the saved boundary wavefield value and
the source (used for concept proofing)
"""
function pressure_reconstruct_forward(path_bnd::Ts, src::Source,
         params::ModelParams) where {Ts <: String}

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
    add_source!(wfd1, src, 1)

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

        # add source to wavefield
        if (params.order+1 <= src.isz <= params.nz-params.order &&
            params.order+1 <= src.isx <= params.nx-params.order)
           add_source!(wfd2, src, it)
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

"""
   Forward reconstruct acoustic pressure field using the saved boundary wavefield value and
the source (used for concept proofing)
"""
function pressure_reconstruct_forward(path_bnd::Ts, srcs::Vector{Source},
         params::ModelParams) where {Ts <: String}

    # length of one-step pressure field
    N = params.nz * params.nx

    # number of sources
    ns= length(srcs)

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
    add_multi_sources!(wfd1, srcs, 1)

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

# """
#    compute recordings with a single source
# """
# function multi_step_forward!(rec::Recordings, src::Source, params::ModelParams)
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
# function multi_step_forward!(rec::Recordings, srcs::Vector{Source}, params::ModelParams)
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
#          bnd::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, params::ModelParams) where {Ti<:Int64, Tv<:AbstractFloat}
#
#     # update vz and vx
#     A_mul_b!(tmp1, rfds.MvzBp, wfd1.p)
#     A_mul_b!(tmp2, rfds.MvxBp, wfd1.p)
#     addition!(wfd2.vz, wfd1.vz, tmp1)
#     addition!(wfd2.vx, wfd1.vx, tmp2)
#
#     # correct for boundary part
#     for i = 1 : length(params.bnd2wfd)
#         j = params.bnd2wfd[i]
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
#     for i = 1 : length(params.bnd2wfd)
#         j = params.bnd2wfd[i]
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
#          params::ModelParams; save_type="pressure")
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
#          ofds::ObsorbFDStencil, params::ModelParams)
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
#          ofds::ObsorbFDStencil, params::ModelParams)
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
#          src::Source, ofds::ObsorbFDStencil, params::ModelParams)
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
#          src::Source, params::ModelParams)
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
