module SeisAcoustic

    # the dependency of this module
    using LinearAlgebra,
          Printf,          # formated print
          DelimitedFiles,  # delimited text file
          SparseArrays,     # finite-difference sparse matrix
          DSP              # smoothing kernal

    # overloading Base function
    import Base.convert, Base.show

    include("dataio/dataio.jl")          # read and write segy data, internally defined regular sampled data (borrow from rsf)
    include("finitediff/finitediff.jl")  # finite-difference method for acoustic wave equation

end # module
