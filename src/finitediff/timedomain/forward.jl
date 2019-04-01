function multiplication!(c::Vector{Tv}, a::Vector{Tv}, b::Vector{Tv}) where {Tv<:AbstractFloat}

    for i = 1 : length(a)
        @inbounds c[i] = a[i] * b[i]
    end
    return nothing
end

function multiplication!(b::Vector{Tv}, a::Vector{Tv}) where {Tv<:AbstractFloat}

    for i = 1 : length(a)
        @inbounds b[i] = b[i] * a[i]
    end
    return nothing
end

function addition!(c::Vector{Tv}, a::Vector{Tv}, b::Vector{Tv}) where {Tv<:AbstractFloat}

    for i = 1 : length(a)
        @inbounds c[i] = a[i]+b[i]
    end
    return nothing
end

function addition!(a::Vector{Tv}, b::Vector{Tv}) where {Tv<:AbstractFloat}
    for i = 1 : length(a)
        @inbounds a[i] = a[i]+b[i]
    end
    return nothing
end

function minus!(a::Vector{Tv}, b::Vector{Tv}, c::Vector{Tv}) where {Tv<:AbstractFloat}
    for i = 1 : length(a)
        @inbounds a[i] = b[i] - c[i]
    end
    return nothing
end

"""
   one forward time-stepping of finite-difference
(rarely used as it take a large amount of memory allocations)
"""
function one_step_forward!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil)

    spt2.vz = ofds.MvzBvz .* spt1.vz + ofds.MvzBp  * (spt1.pz+spt1.px)
    spt2.vx = ofds.MvxBvx .* spt1.vx + ofds.MvxBp  * (spt1.pz+spt1.px)
    spt2.pz = ofds.MpzBpz .* spt1.pz + ofds.MpzBvz *  spt2.vz
    spt2.px = ofds.MpxBpx .* spt1.px + ofds.MpxBvx *  spt2.vx

    return nothing
end

"""
   one forward time-stepping of finite-difference
"""
function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, ofds::ObsorbFDStencil{Ti,Tv},
         tmp1::Vector{Tv}, tmp2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    addition!(tmp1, spt1.pz, spt1.px)               # p1 = pz1+px1

    # update vz
    multiplication!(spt2.vz, ofds.MvzBvz, spt1.vz)  # vz2 = MvzBvz * vz1
    A_mul_b!(tmp2, ofds.MvzBp, tmp1)                # tmp2= MvzBp  * p1
    addition!(spt2.vz, tmp2)                        # vz2 = vz2    + tmp2

    # update vx
    multiplication!(spt2.vx, ofds.MvxBvx, spt1.vx)  # vx2 = MvxBvx * vx1
    A_mul_b!(tmp2, ofds.MvxBp, tmp1)                # tmp2= MvxBp  * p1
    addition!(spt2.vx, tmp2)                        # vx2 = vx2    + tmp2

    # update pz
    multiplication!(spt2.pz, ofds.MpzBpz, spt1.pz)  # pz2 = MpzBpz * pz1
    A_mul_b!(tmp1, ofds.MpzBvz, spt2.vz)            # tmp1= MpzBvz * vz2
    addition!(spt2.pz, tmp1)                        # pz2 = pz2    + tmp1

    # update px
    multiplication!(spt2.px, ofds.MpxBpx, spt1.px)  # pzx = MpxBpx * px1
    A_mul_b!(tmp1, ofds.MpxBvx, spt2.vx)            # tmp1= MpxBvx * vx2
    addition!(spt2.px, tmp1)                        # px2 = px2    + tmp1

    return nothing
end

"""
   one step Forward reconstruction of wave field from the boundary values
"""
function one_step_forward!(wfd2::Wavefield, wfd1::Wavefield, rfds::RigidFDStencil, it::Ti,
         bnd::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, params::ModelParams) where {Ti<:Int64, Tv<:AbstractFloat}

    # update vz and vx
    A_mul_b!(tmp1, rfds.MvzBp, wfd1.p)
    A_mul_b!(tmp2, rfds.MvxBp, wfd1.p)
    addition!(wfd2.vz, wfd1.vz, tmp1)
    addition!(wfd2.vx, wfd1.vx, tmp2)

    # correct for boundary part
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd2.vz[j] = bnd.vz[i,it]
        wfd2.vx[j] = bnd.vx[i,it]
    end

    # update pressure
    A_mul_b!(tmp1, rfds.MpxBvx, wfd2.vx)
    A_mul_b!(tmp2, rfds.MpzBvz, wfd2.vz)
    addition!(tmp1, tmp2)
    addition!(wfd2.p, wfd1.p, tmp1)

    # correct for boundary part
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd2.p[j] = bnd.p[i,it]
    end

    return nothing
end

"""
   one step backward reconstruction of wave field from the boundary values
"""
function one_step_backward!(wfd1::Wavefield, wfd2::Wavefield, rfds::RigidFDStencil, it::Ti,
         bnd::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, params::ModelParams) where {Ti<:Int64, Tv<:AbstractFloat}

    # wfd1.p  = wfd2.p  - ofds.MpxBvx * wfd2.vx - ofds.MpzBvz * wfd2.vz
    # wfd1.p[bound.index]  = bound.p[:, it]
    A_mul_b!(tmp1, rfds.MpxBvx, wfd2.vx)
    A_mul_b!(tmp2, rfds.MpzBvz, wfd2.vz)
    addition!(tmp1, tmp2)
    minus!(wfd1.p, wfd2.p, tmp1)

    # correct for boundary part of pressure
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd1.p[j] = bnd.p[i,it]
    end

    # wfd1.vx = wfd2.vx - ofds.MvxBp  * wfd1.p
    # wfd1.vz = wfd2.vz - ofds.MvzBp  * wfd1.p
    # wfd1.vx[bound.index] = bound.vx[:, it]
    # wfd1.vz[bound.index] = bound.vz[:, it]
    A_mul_b!(tmp1, ofds.MvxBp, wfd1.p)
    A_mul_b!(tmp2, ofds.MvzBp, wfd1.p)
    minus!(wfd1.vx, wfd2.vx, tmp1)
    minus!(wfd1.vz, wfd2.vz, tmp2)

    # correct for boundary part of particle velocity
    for i = 1 : length(params.bnd2wfd)
        j = params.bnd2wfd[i]
        wfd1.vz[j] = bnd.vz[i,it]
        wfd1.vx[j] = bnd.vx[i,it]
    end

    return nothing
end

"""
save snapshot or full wave field or pressure field to hard drive, inject single source
save_type="snapshot" : vz, vx, pz, px include PML part,
save_type="wavefield": vz, vx, pz+px  without PML part,
save_type="pressure" : p=pz+px without PML bounary part.
"""
function multi_step_forward(path::String, src::Source, ofds::ObsorbFDStencil,
         params::ModelParams; save_type="pressure")

    # initialize some variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # save snapshot
    if save_type == "snapshot"
       hdr = snapshot_header(par)
       fid = write_RSheader(path, hdr)
       append_one_snapshot(fid, spt1)

    # save wavefield
    elseif save_type == "wavefield"
       hdr = wavefield_header(par)
       fid = write_RSheader(path, hdr)
       append_one_wavefield(fid, spt1, params)

    # save pressure
    elseif save_type == "pressure"
       hdr = pressure_header(par)
       fid = write_RSheader(path, hdr)
       append_one_pressure(fid, spt1, params)
    else
       error("save_type only support snapshot, wavefield or pressure")
    end

    # loop over time-stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)
        copy_snapshot!(spt1, spt2)

        if save_type == "snapshot"
           append_one_snapshot(fid, spt1)
        elseif save_type == "wavefield"
           append_one_wavefield(fid, spt1, params)
        elseif save_type == "pressure"
           append_one_pressure(fid, spt1, params)
        end
    end

    close(fid)
    return nothing
end

"""
   compute wavefield boundary at each time step and the wavefield at last time step
save them to hard drive.
"""
function get_boundary_wavefield(path_wfd::String, path_bnd::String,
         src::Source, ofds::ObsorbFDStencil, params::ModelParams)

    # one extra time step
    nt = par.nt + 1

    # initialize variables for time stepping
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # save the boundary of the first wavefield
    hdr = boundary_header(params)
    fid = write_RSheader(path_bnd, hdr)
    append_one_boundary(fid, spt1, params)

    # loop over time stepping
    for it = 2 : nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)

        # save boundary of wavefield
        append_one_boundary(fid, spt2, params)

        # prepare for next iteration
        copy_snapshot!(spt1, spt2)
    end

    # finish saving all the boundary values, close the file
    close(fid)

    # save the wavefield at the final time step
    wfd = sample_spt2wfd(spt1, params)
    write_wavefield(path_wfd, wfd)

    return nothing
end

"""
   compute recordings with a single source
"""
function multi_step_forward!(rec::Recordings, src::Source,
         ofds::ObsorbFDStencil, params::ModelParams)

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # obtain the first time sample of recordings
    sample_spt2rec!(rec, spt1, 1)

    # loop over time stepping
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)
        sample_spt2rec!(rec, spt2, it)

        # prepare for next time step
        copy_snapshot!(spt1, spt2)
    end

    return nothing
end

"""
   compute recordings with multiple simultaneous source, the only difference is that the input
is a vector of source
"""
function multi_step_forward!(rec::Recordings, srcs::Vector{Source},
         ofds::ObsorbFDStencil, params::ModelParams)

    # initialize intermediate variables
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp1 = zeros(params.data_format, params.Nz * params.Nx)
    tmp2 = zeros(params.data_format, params.Nz * params.Nx)

    # add multi-sources to the first snapshot
    add_multi_sources!(spt1, srcs, 1)

    # sample the snapshot to fill recordings
    sample_spt2rec!(rec, spt1, 1)

    # loop over time
    for it = 2 : params.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_multi_sources!(spt2, srcs, it)
        sample_spt2rec!(rec, spt2, it)

        # prepare for the next time step
        copy_snapshot!(spt1, spt2)
    end

    return nothing
end

"""
   Forward reconstruct acoustic pressure field using the saved boundary wavefield value and
the source (used for concept proofing)
"""
function pressurefield_reconstruct_forward(bnd::WavefieldBound, rfds::RigidFDStencil,
         src::Source, params::ModelParams)

    # length of one-step pressure field
    N = params.nz * params.nx

    # decide the number of boundary layers
    order = length(params.fd_coefficients)

    # initialize intermediate variables
    wfd1 = Wavefield(params)
    wfd2 = Wavefield(params)
    tmp1 = zeros(params.data_format, N)
    tmp2 = zeros(params.data_format, N)

    # allocate memory for saving pressure field
    pre = zeros(par.data_format, N * params.nt)

    # add source to the first wavefield
    add_source!(wfd1, src, 1)

    # save pressure field
    copyto!(pre, 1, wfd1.p, 1, N)

    # the start index for saving the next pressure field
    idx_o = N+1

    # loop over time stepping
    for it = 2 : params.nt

        # forward time steping and correcting boundaries
        one_step_forward!(wfd2, wfd1, rfds, it, bnd, tmp1, tmp2, params)

        # add source to wavefield
        if order+1 <= src.isz <= par.nz-order && order+1 <= src.isx <= par.nx-order
           add_source!(wfd2, src, it)
        end

        # save the pressure field
        copyto!(pre, idx_o, wfd2.p, 1, N)
        idx_o = idx_o + 1

        # prepare for next step
        copy_wavefield!(wfd1, wfd2)
    end

    return reshape(pre, params.nz, params.nx, params.nt)
end
