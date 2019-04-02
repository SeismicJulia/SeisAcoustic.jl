module SeisAcoustic

    # the dependency of this module
    using LinearAlgebra,
          Printf,           # formated print
          DelimitedFiles,   # delimited text file
          SparseArrays,     # finite-difference sparse matrix
          DSP,              # smoothing kernal
          FFTW,             # fourier transform
          BinDeps           # manage binary dependencies

    # the absolute path to the binary dependency
    const spmatveclib  = abspath(joinpath(splitdir(Base.source_path())[1],"..","deps","builds","spmatvec.so"))

    # overloading Base function
    import Base.convert, Base.show

    include("dataio/dataio.jl")          # read and write segy data, internally defined regular sampled data (borrow from rsf)
    include("finitediff/finitediff.jl")  # finite-difference method for acoustic wave equation

end # module
