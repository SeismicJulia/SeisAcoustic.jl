# ==================== struct of model parameters ===============
export get_gamma_profile,
       Mvz,
       Rvz,
       Mvx,
       Rvx,
       Mpz,
       Rpz,
       Mpx,
       Rpx,
       TdParams,
       show,
       model_padding,
       model_smooth,
       density2buoyancy,
       vp2lambda,
       finite_difference_coefficients,
# ==================== Snapshot =================
       Snapshot,
       Wavefield,
       sample_spt2wfd,
       sample_adjoint2pre!,
       snapshot_header,
       append_one_snapshot,
       read_one_snapshot,
# ==================== Wavefield =================
       wavefield_header,
       append_one_wavefield,
       write_wavefield,
       read_one_wavefield,
# ==================== pressure =================
       pressure_header,
       append_one_pressure,
       read_one_pressure,
       copy_snapshot!,
       copy_wavefield!,
# ==================== Boundary =================
       WavefieldBound,
       boundary_header,
       write_boundary,
       append_one_boundary,
       read_boundary,
       read_one_boundary!,
       l2norm_snapshot,
       minus_snapshot,
       add_snapshot,
       reverse_order,
# # ==================== Source =================
       Source,
       ricker,
       convert2miniphase,
       tapered_sinc,
       scaled_sinc,
       get_multi_sources,
       add_source!,
       subtract_source!,
       time_range_source,
       write_source,
       read_source,
       source_isequal,
# # ==================== Recordings ===============
       Recordings,
       write_recordings,
       read_recordings,
       sample_spt2rec!,
       inject_rec2spt!,
       recordings_isequal,
# ============= forward time stepping =============
       A_mul_b!,
       At_mul_b!,
# ============= forward time stepping =============
       vertical_partial_derivative!,
       horizontal_partial_derivative!,
       one_step_forward!,
       multi_step_forward,
       pressure_reconstruct_forward,
# ============= adjoint time stepping =============
       one_step_adjoint!,
       multi_step_adjoint!,
       one_step_backward!,
       pressure_reconstruct_backward

include("pml.jl")
include("fdstencil.jl")
include("tdparams.jl")
include("snapshot.jl")
include("source.jl")
include("recordings.jl")
include("spmatvec.jl")
include("forward.jl")
include("adjoint.jl")



# # ==================== struct of model parameters ===============
# export ModelParams,
#        show,
#        model_padding,
#        model_smooth,
#        finite_difference_coefficients,
# # ==================== PML ======================
#        PMLCoefficients,
# # ==================== FD stencil ===============
#        ObsorbFDStencil,
#        RigidFDStencil,
# # ==================== Snapshot =================
#        Snapshot,
#        Wavefield,
#        sample_spt2wfd,
#        snapshot_header,
#        append_one_snapshot,
#        read_one_snapshot,
# # ==================== Wavefield =================
#        wavefield_header,
#        append_one_wavefield,
#        write_wavefield,
#        read_one_wavefield,
# # ==================== pressure =================
#        pressure_header,
#        append_one_pressure,
#        read_one_pressure,
#        copy_snapshot!,
#        copy_wavefield!,
# # ==================== Boundary =================
#        WavefieldBound,
#        boundary_header,
#        write_boundary,
#        append_one_boundary,
#        read_boundary,
#        read_one_boundary,
#        l2norm_snapshot,
#        minus_snapshot,
#        add_snapshot,
# # # ==================== Source =================
#        Source,
#        ricker,
#        convert2miniphase,
#        tapered_sinc,
#        scaled_sinc,
#        get_multi_sources,
#        add_source!,
#        add_multi_sources!,
#        subtract_source!,
#        time_range_multisources,
# # # ==================== Recordings ===============
#        Recordings,
#        write_recordings,
#        read_recordings,
#        sample_spt2rec!,
#        inject_rec2spt!,
# # ============= forward time stepping =============
#        A_mul_b!,
#        At_mul_b!,
# # ============= forward time stepping =============
#        one_step_forward!,
#        one_step_backward!,
#        multi_step_forward,
#        multi_step_forward!,
#        get_boundary_wavefield,
#        pressure_reconstruct_forward,
# # ============= adjoint time stepping =============
#        one_step_adjoint!,
#        one_step_backward,
#        multi_step_adjoint,
#        pressure_reconstruct_backward
#
#
# include("modelparams.jl")
# include("pml.jl")
# include("fdstencil.jl")
# include("snapshot.jl")
# include("source.jl")
# include("recordings.jl")
# include("spmatvec.jl")
# include("forward.jl")
# include("adjoint.jl")
