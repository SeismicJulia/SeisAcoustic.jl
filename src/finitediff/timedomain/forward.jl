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
   one forward time-stepping of finite-difference (rarely used)
"""
function one_step_forward!(spt2::Snapshot, spt1::Snapshot, ofds::ObsorbFDStencil)

    spt2.vz = ofds.MvzBvz .* spt1.vz + ofds.MvzBp  * (spt1.pz+spt1.px)
    spt2.vx = ofds.MvxBvx .* spt1.vx + ofds.MvxBp  * (spt1.pz+spt1.px)
    spt2.pz = ofds.MpzBpz .* spt1.pz + ofds.MpzBvz *  spt2.vz
    spt2.px = ofds.MpxBpx .* spt1.px + ofds.MpxBvx *  spt2.vx

    return nothing
end

"""
   in-place one forward time-stepping of finite-difference
"""
function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, ofds::ObsorbFDStencil{Ti,Tv},
         tmp1::Vector{Tv}, tmp2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    addition!(tmp1, spt1.pz, spt1.px)               # p = pz+px

    # update vz
    multiplication!(spt2.vz, ofds.MvzBvz, spt1.vz)  # vz2 = MvzBvz * vz
    AmulB!(tmp2, ofds.MvzBp, tmp1)                  # MvzBp * (pz+px)
    addition!(spt2.vz, tmp2)

    # update vx
    multiplication!(spt2.vx, ofds.MvxBvx, spt1.vx)
    AmulB!(tmp2, ofds.MvxBp, tmp1)
    addition!(spt2.vx, tmp2)

    # update pz
    multiplication!(spt2.pz, ofds.MpzBpz, spt1.pz)
    AmulB!(tmp1, ofds.MpzBvz, spt2.vz)
    addition!(spt2.pz, tmp1)

    # update px
    multiplication!(spt2.px, ofds.MpxBpx, spt1.px)
    AmulB!(tmp1, ofds.MpxBvx, spt2.vx)
    addition!(spt2.px, tmp1)

    return nothing
end

"""
   one step Forward reconstruction of wave field from the boundary values
"""
function one_step_forward!(wfd2::Wavefield, wfd1::Wavefield, ofds::RigidFiniteDiffMatrix, it::Ti,
         bound::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, par::PhysicalModel) where {Ti<:Int64, Tv<:AbstractFloat}

    # update vz and vx
    AmulB!(tmp1, ofds.MvzBp, wfd1.p)
    AmulB!(tmp2, ofds.MvxBp, wfd1.p)
    addition!(wfd2.vz, wfd1.vz, tmp1)
    addition!(wfd2.vx, wfd1.vx, tmp2)

    # correct for boundary part
    for i = 1 : length(par.index3)
        j = par.index3[i]
        wfd2.vz[j] = bound.vz[i,it]
        wfd2.vx[j] = bound.vx[i,it]
    end

    # update pressure
    AmulB!(tmp1, ofds.MpxBvx, wfd2.vx)
    AmulB!(tmp2, ofds.MpzBvz, wfd2.vz)
    addition!(tmp1, tmp2)
    addition!(wfd2.p, wfd1.p, tmp1)

    # correct for boundary part
    for i = 1 : length(par.index3)
        j = par.index3[i]
        wfd2.p[j] = bound.p[i,it]
    end

    return nothing
end

"""
   one step backward reconstruction of wave field from the boundary values
"""
function one_step_backward!(wfd1::Wavefield, wfd2::Wavefield, ofds::RigidFiniteDiffMatrix, it::Ti,
         bound::WavefieldBound, tmp1::Vector{Tv}, tmp2::Vector{Tv}, par::PhysicalModel) where {Ti<:Int64, Tv<:AbstractFloat}

    # wfd1.p  = wfd2.p  - ofds.MpxBvx * wfd2.vx - ofds.MpzBvz * wfd2.vz
    # wfd1.p[bound.index]  = bound.p[:, it]
    AmulB!(tmp1, ofds.MpxBvx, wfd2.vx)
    AmulB!(tmp2, ofds.MpzBvz, wfd2.vz)
    addition!(tmp1, tmp2)
    minus!(wfd1.p, wfd2.p, tmp1)

    for i = 1 : length(par.index3)
        j = par.index3[i]
        wfd1.p[j] = bound.p[i,it]
    end

    # wfd1.vx = wfd2.vx - ofds.MvxBp  * wfd1.p
    # wfd1.vz = wfd2.vz - ofds.MvzBp  * wfd1.p
    # wfd1.vx[bound.index] = bound.vx[:, it]
    # wfd1.vz[bound.index] = bound.vz[:, it]
    AmulB!(tmp1, ofds.MvxBp, wfd1.p)
    AmulB!(tmp2, ofds.MvzBp, wfd1.p)
    minus!(wfd1.vx, wfd2.vx, tmp1)
    minus!(wfd1.vz, wfd2.vz, tmp2)

    for i = 1 : length(par.index3)
        j = par.index3[i]
        wfd1.vz[j] = bound.vz[i,it]
        wfd1.vx[j] = bound.vx[i,it]
    end

    return nothing

end

"""
save snapshot // full wave field // stress field to hard drive, inject single source
save_type="snapshot" : vz, vx, pz, px include PML part,
save_type="wavefield": vz, vx, pz+px  without PML part,
save_type="pressyre" : p=pz+px without PML bounary part.
"""
function multi_step_forward(path::String, src::Source, ofds::FiniteDiffMatrix,
         par::PhysicalModel; save_type="snapshot")

    # initialize some variables
    spt1 = Snapshot(par)
    spt2 = Snapshot(par)
    tmp1 = zeros(spt1.vz)
    tmp2 = zeros(spt1.vx)

    if src.ot == 0.0
       add_source!(spt1, src, 1)
    end

    if save_type == "snapshot"
       hdr = snapshot_header(par)
       fid = write_USdata(path, hdr)
       append_snapshot(fid, spt1)
    elseif save_type == "wavefield"
       hdr = wavefield_header(par)
       fid = write_USdata(path, hdr)
       tmp = zeros(par.data_format, par.nz*par.nx)
       append_wavefield(fid, spt1, par, tmp)
    elseif save_type == "pressure"
       hdr = pressure_header(par)
       fid = write_USdata(path, hdr)
       tmp = zeros(par.data_format, par.nz*par.nx)
       append_pressure(fid, spt1, par, tmp)
    end

    for it = 2 : par.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)
        copy_snapshot!(spt1, spt2)

        if save_type == "snapshot"
           append_snapshot(fid, spt1)
        elseif save_type == "wavefield"
           append_wavefield(fid, spt1, par, tmp)
        elseif save_type == "pressure"
           append_pressure(fid, spt1, par, tmp)
        end
    end
    close(fid)
end

"""
   output pressure wavefield
"""
function multi_step_forward(src::Source, ofds::FiniteDiffMatrix, par::PhysicalModel)

    # initialize some variables
    spt1 = Snapshot(par)
    spt2 = Snapshot(par)
    tmp1 = zeros(spt1.vz)
    tmp2 = zeros(spt1.vx)

    # the output
    N = par.nz * par.nx
    wfd  = zeros(par.data_format, N, par.nt)

    if src.ot == 0.0
       add_source!(spt1, src, 1)
    end

    function spt2wfd!(wfd::Matrix{Tv}, spt::Snapshot, it::Ti,
             N::Ti, par::PhysicalModel) where{Tv<:AbstractFloat, Ti<:Int64}

        for i = 1 : N
            wfd[i,it] = spt.pz[par.index1[i]] + spt.px[par.index1[i]]
        end
    end
    spt2wfd!(wfd, spt1, 1, N, par)

    for it = 2 : par.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)
        copy_snapshot!(spt1, spt2)
        spt2wfd!(wfd, spt1, it, N, par)
    end

    return wfd
end

function save_boundary!(bnd::WavefieldBound, spt::Snapshot, index::Vector{Ti}, it::Ti) where {Ti<:Int64}

    for i = 1 : length(index)
        j = index[i]
        bnd.vz[i,it] = spt.vz[j]
        bnd.vx[i,it] = spt.vx[j]
        bnd.p[i ,it] = spt.pz[j] + spt.px[j]
    end
    return nothing
end

"""
   compute the boundary of each wavefield and the last wavefield
"""
function get_wavefield_bound(src::Source, ofds::FiniteDiffMatrix, par::PhysicalModel;
         save_flag=false, path1="NULL", path2="NULL")

    nt= par.nt + 1
    bnd = WavefieldBound(par)

    # initialize
    spt1 = Snapshot(par)
    spt2 = Snapshot(par)
    tmp1 = zeros(spt1.vz)
    tmp2 = zeros(spt1.vx)
    one_over_two = par.data_format(1./2.0)
    one          = par.data_format(1.    )

    # initialize the array to save result
    N = par.nz * par.nx
    src_strength = zeros(par.data_format, N)
    mp   = zeros(par.data_format, N)

    # add source
    if src.ot == 0.0
       add_source!(spt1, src, 1)
    end
    save_boundary!(bnd, spt1, par.index2, 1)

    # compute the strength of source side wavefield (the time direvative of wavefield)
    function compute_source_strength!(src_strength::Vector{Tv}, N::Ti, spt1::Snapshot,
             spt2::Snapshot, mp::Vector{Tv}, one_over_two::Tv) where {Ti<:Int64, Tv<:AbstractFloat}

        for i = 1 : N
            j = par.index1[i]
            tmp = (spt2.pz[j]+spt2.px[j]-mp[i]) * one_over_two
            src_strength[i] = src_strength[i] + tmp * tmp
            mp[i] = spt1.pz[j] + spt1.px[j]
        end
    end

    for it = 2 : nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)

        # save boundary of wavefield
        save_boundary!(bnd, spt2, par.index2, it)

        # compute the strenght of source-side wavefield
        compute_source_strength!(src_strength, N, spt1, spt2, mp, one_over_two)

        # prepare for next iteration
        copy_snapshot!(spt1, spt2)
    end

    # save the wavefield at the final time step
    wfd = sample_spt2wavefield(spt1, par)

    if save_flag
       if path1 == "NULL" || path2 == "NULL"
          error("does not provide direcotry to save boundary wavefield or source strength")
       end
       write_bound_wavefield(path1, bnd, wfd)

       hdr = UniformSampleHeader(src_strength, title="source-side wavefield strength")
       write_USdata(path2, hdr, src_strength)
       return nothing

    else
       return bnd, wfd, src_strength

    end
end

# output record or OBC record, inject single source
function multi_step_forward(irz::Vector{Ti}, irx::Vector{Ti}, src::Source, ofds::FiniteDiffMatrix,
         par::PhysicalModel; record_type="record") where {Ti<:Int64}

    # initialize output
    if record_type == "record"
       rec = Record(src.isz, src.isx, irz, irx, par)
    elseif record_type == "recordOBC"
       rec = RecordOBC(src.isz, src.isx, irz, irx, par)
    end

    # initialize intermediate variables
    spt1 = Snapshot(par)
    spt2 = Snapshot(par)
    tmp1 = zeros(spt1.vz)
    tmp2 = zeros(spt1.vx)

    # the first snapshot
    if src.ot == 0.0
       add_source!(spt1, src, 1)
    end

    if record_type == "record"
       sample_spt2record!(rec, spt1, 1)
    elseif record_type == "recordOBC"
       sample_spt2recordOBC!(rec, spt1, 1)
    end

    for it = 2 : par.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_source!(spt2, src, it)
        copy_snapshot!(spt1, spt2)

        if record_type == "record"
           sample_spt2record!(rec, spt1, it)
        elseif record_type == "recordOBC"
           sample_spt2recordOBC!(rec, spt1, it)
        end

    end

    return rec
end

# output record or OBC record, inject multi-sources
function multi_step_forward(irz::Vector{Ti}, irx::Vector{Ti}, srcs::Vector{Source}, ofds::FiniteDiffMatrix,
         par::PhysicalModel; record_type="record") where {Ti<:Int64}

    # initialize output
    if record_type == "record"
       rec = Record(0, 0, irz, irx, par)
    elseif record_type == "recordOBC"
       rec = RecordOBC(0, 0, irz, irx, par)
    end

    # initialize intermediate variables
    spt1 = Snapshot(par)
    spt2 = Snapshot(par)
    tmp1 = zeros(spt1.vz)
    tmp2 = zeros(spt1.vx)

    # the first snapshot
    (tl, tu) = time_range_multisources(srcs)
    if tl == 0.0
       add_multi_sources!(spt1, srcs, 1)
    end

    if record_type == "record"
       sample_spt2record!(rec, spt1, 1)
    elseif record_type == "recordOBC"
       sample_spt2recordOBC!(rec, spt1, 1)
    end

    for it = 2 : par.nt

        one_step_forward!(spt2, spt1, ofds, tmp1, tmp2)
        add_multi_sources!(spt2, srcs, it)
        copy_snapshot!(spt1, spt2)

        if record_type == "record"
           sample_spt2record!(rec, spt1, it)
        elseif record_type == "recordOBC"
           sample_spt2recordOBC!(rec, spt1, it)
        end

    end

    return rec
end

"""
   Reconstruct acoustic-wavefield forward using the boundary wavefield value and
the source, it equvialent to forward modeling without PML boundary condition.
"""
function wavefield_reconstruct_forward(bnd::WavefieldBound, ofdsT::RigidFiniteDiffMatrix,
         src::Source, par::PhysicalModel)

    N = par.nz * par.nx
    order = length(par.fd)
    wfd1 = Wavefield(par)
    wfd2 = Wavefield(par)
    tmp1 = zeros(par.data_format, N)
    tmp2 = zeros(par.data_format, N)

    add_source!(wfd1, src, 1)

    # allocate memory for saving the wave
    pre = zeros(par.data_format, N, par.nt)
    pre[:,1] = wfd1.p

    for it = 2 : par.nt
        one_step_forward!(wfd2, wfd1, ofdsT, it, bnd, tmp1, tmp2, par)
        if order+1 <= src.isz <= par.nz-order && order+1 <= src.isx <= par.nx-order
           add_source!(wfd2, src, it)
        end
        pre[:,it] = wfd2.p
        copy_wavefield!(wfd1, wfd2)
    end
    return pre
end


# return shot or shot3, inject stress which read from slow memory
# function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String , ofds::ofds; tmax=1.0, otype="p")
#     nz = ofds.nz ; nx = ofds.nx; ext= ofds.ext; iflag = ofds.iflag; dt = ofds.dt;
#     nt = round(Int64, tmax/dt)+1
#     (nz1, nx1, ext1, iflag1, dt1, nt1) = InfoStress(path)
#     if otype == "p"
#        shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     elseif otype == "vp"
#        shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     end
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     addStress2Spt!(spt1, path)
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp)
#     if otype == "p"
#        spt2shot!(shot, spt1)
#     elseif otype == "vp"
#        spt2shot3!(shot, spt1)
#     end
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         if it <= nt1
#            addStress2Spt!(spt2, path)
#         end
#         copySnapShot!(spt1, spt2)
#         if otype == "p"
#            spt2shot!(shot, spt1)
#         elseif otype == "vp"
#            spt2shot3!(shot, spt1)
#         end
#     end
#     return shot
# end
#
# # return shot or shot3, do born approximation
# function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, I::Array{Float64,2}, path::String , ofds::ofds; tmax=1.0, otype="p")
#     nz = ofds.nz ; nx = ofds.nx; ext= ofds.ext; iflag = ofds.iflag; dt = ofds.dt;
#     nt = round(Int64, tmax/dt) + 1
#     (nz1, nx1, ext1, iflag1, dt1, nt1) = InfoStress(path)
#     if otype == "p"
#        shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     elseif otype == "vp"
#        shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     end
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     addReflection2Spt!(spt1, I, path)
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp);
#     if otype == "p"
#        spt2shot!(shot, spt1)
#     elseif otype == "vp"
#        spt2shot3!(shot, spt1)
#     end
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         if it <= nt1
#            addReflection2Spt!(spt2, I, path)
#         end
#         copySnapShot!(spt1, spt2)
#         if otype == "p"
#            spt2shot!(shot, spt1)
#         elseif otype == "vp"
#            spt2shot3!(shot, spt1)
#         end
#     end
#     return shot
# end
#
# #  output shot or shot3, input wavelet and the distribution of source(used for source location)
# function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, w::Array{Float64,1}, dis::Array{Float64,2}, ofds::ofds; tmax=1.0, otype="p")
#     nz = ofds.nz ; nx = ofds.nx; ext= ofds.ext; iflag = ofds.iflag; dt = ofds.dt;
#     nt = round(Int64, tmax/dt) + 1
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nx +   ext
#     end
#     Nx  = nx + 2*ext
#     nst = length(w)
#     if otype == "p"
#        shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     elseif otype == "vp"
#        shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     end
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     AddMultiSources!(spt1, w, dis)
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp);
#     if otype == "p"
#        spt2shot!(shot, spt1)
#     elseif otype == "vp"
#        spt2shot3!(shot, spt1)
#     end
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         if it <= nst
#            AddMultiSources!(spt2, w, dis)
#         end
#         copySnapShot!(spt1, spt2)
#         if otype == "p"
#            spt2shot!(shot, spt1)
#         elseif otype == "vp"
#            spt2shot3!(shot, spt1)
#         end
#     end
#     return shot
# end
#
# # save snapshot // full wave field // stress field to slow memory, inject single source
# function MultiStepForward(src::Source, ofds::ofds, interval::Int64; tmax=1.0)
#
#     etype = eltype(ofds.MvzBp.nzval)
#
#     nz = ofds.nz ; nx = ofds.nx;
#     ext= ofds.ext; iflag = ofds.iflag;
#     dt = ofds.dt;
#
#     if iflag == 1
#        Nz = nz + 2*ext
#        zupper = ext
#     elseif iflag == 2
#        Nz = nz +   ext
#        zupper = 0
#     end
#     Nx = nx + 2*ext
#
#     nt = round(Int64, tmax/dt)+1
#     tl = src.ot; tu = tl + (src.nt-1)*dt;
#
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#
#     AddSource!(spt1, src)
#
#     wfd = spt1.pz + spt1.px; wfd = reshape(wfd, Nz, Nx);
#     wfd = wfd[zupper+1:zupper+nz, ext+1:ext+nx]
#
#     nt1 = floor(Int64, tmax/(dt*interval)) + 1;
#     result = zeros(etype, nz, nx, nt1);
#
#     idx = 1;
#     result[:,:,idx] = wfd;
#     idx = idx + 1;
#
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp);
#
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         if tl <= (it-1)*dt <= tu
#            AddSource!(spt2, src)
#         end
#         copySnapShot!(spt1, spt2)
#         if mod(it-1, interval) == 0
#            wfd = spt1.pz + spt1.px; wfd = reshape(wfd, Nz, Nx);
#            wfd = wfd[zupper+1:zupper+nz, ext+1:ext+nx]
#            result[:,:,idx] = wfd; idx = idx + 1;
#         end
#     end
#     return result
# end
#
# # inject encoded sources and the output is the whole wave-field
# function MultiStepForward(srcs::Vector{Source}, ofds::ofds, interval::Int64; tmax=1.0)
#
#     nz = ofds.nz ; nx    = ofds.nx;
#     ext= ofds.ext; iflag = ofds.iflag;
#     dt = ofds.dt;
#
#     if iflag == 1
#        Nz = nz + 2*ext
#        zupper = ext
#     elseif iflag == 2
#        Nz = nz +   ext
#        zupper = 0
#     end
#     Nx = nx + 2*ext
#
#     nt = round(Int64, tmax/dt)+1
#     (tl, tu) = SrcRange(srcs)
#
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#
#     AddMultiSources!(spt1, srcs)
#
#     wfd = spt1.pz + spt1.px; wfd = reshape(wfd, Nz, Nx);
#     wfd = wfd[zupper+1:zupper+nz, ext+1:ext+nx]
#
#     nt1 = floor(Int64, tmax/(dt*interval)) + 1;
#     result = zeros(Float32, nz, nx, nt1);
#
#     idx = 1;
#     result[:,:,idx] = wfd;
#     idx = idx + 1;
#
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp);
#
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         if tl <= (it-1)*dt <= tu
#            AddMultiSources!(spt2, srcs)
#         end
#         copySnapShot!(spt1, spt2)
#         if mod(it-1, interval) == 0
#            wfd = spt1.pz + spt1.px; wfd = reshape(wfd, Nz, Nx);
#            wfd = wfd[zupper+1:zupper+nz, ext+1:ext+nx]
#            result[:,:,idx] = wfd; idx = idx + 1;
#         end
#     end
#     return result
# end

# output one shot, inject windowed source cube
# function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, wsc::WSrcCube, ofds::ofds; tmax=1.0)
#     nz = ofds.nz ; nx = ofds.nx; ext= ofds.ext; iflag = ofds.iflag; dt = ofds.dt;
#     nt = round(Int64, tmax/dt)+1
#     shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     AddWsc2spt!(spt1, wsc)
#     spt2shot!(shot, spt1)
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp);
#     for it = 2 : wsc.nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         AddWsc2spt!(spt2, wsc)
#         copySnapShot!(spt1, spt2)
#         spt2shot!(shot, spt1)
#     end
#     for it = wsc.nt+1 : nt
#         OneStepForward!(spt2, spt1, ofds, tmp, tmp1)
#         copySnapShot!(spt1, spt2)
#         spt2shot!(shot, spt1)
#     end
#     return shot
# end
