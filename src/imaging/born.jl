"""
   Add the virtual source to pressure (just add to pz), the source-side wavefield
of all time-steps are saved in memory.
"""
function add_virtual_source!(spt::Snapshot, m::Vector{Tv}, wfd1::Wavefield, wfd2::Wavefield,
                             params::TdParams) where {Tv<:AbstractFloat}

    # total number of element
    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]
        dpdt = 2.0 * (wfd2.p[i] - wfd1.p[i]) / params.vel[i]
        spt.px[j] = spt.px[j] +  dpdt * m[i]
    end

    return nothing
end

"""
   the first sample of recordings can be anything, the scattered wavefield is computed
from the second time step.
"""
function born_approximation_forward!(rec::Recordings, m::Vector{Tv}, path_bnd::Ts,
                                     src::Source, params::TdParams) where {Ts<:String, Tv<:AbstractFloat}

    # model length
    N = params.nz * params.nx

    # allocate memory for computing adjoint wavefield
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # allocate memory for reconstructing source-side wavefield backward
    wfd1 = Wavefield(params); add_source!(wfd1, src, 1); # get the wavefield at first time step
    wfd2 = Wavefield(params)
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # forward time stepping
    for it = 2 : params.nt

        # time stepping of scattered wavefield
        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # reconstruct source-side wavefield
        read_one_boundary!(bnd, fid_bnd, it, params)
        one_step_forward!(wfd2, wfd1, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # subtract source term if it being added during update boundary part
        if (src.isz <= params.order || src.isz > params.nz-params.order ||
            src.isx <= params.order || src.isx > params.nx-params.order)

            subtract_source!(wfd2, src, it)
        end

        # add the virtual source
        add_virtual_source!(spt2, m, wfd1, wfd2, params)

        # sampling scatter wavefield
        sample_spt2rec!(rec, spt2, it)

        # prepare for the next step backward reconstruction
        add_source!(wfd2, src, it)
        copy_wavefield!(wfd1, wfd2)
        copy_snapshot!(spt1, spt2)

    end

    # close the boundary value file
    close(fid_bnd)

    return nothing
end
