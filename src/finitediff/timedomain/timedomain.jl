# ==================== struct of model parameters ===============
export ModelParams,
       show,
       model_padding,
       model_smooth,
       finite_difference_coefficients,
# ==================== PML ======================
       PMLCoefficients,
# ==================== FD stencil ===============
       ObsorbFDStencil,
       RigidFDStencil


include("modelparams.jl")
include("pml.jl")
include("fdstencil.jl")
