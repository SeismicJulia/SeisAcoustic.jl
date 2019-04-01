"""
   y = A * x, A is a sparse matrix in CSC format, x is a vector, y is updated in-place
can be used for both Float32 or Float64.
"""
function A_mul_b!(y::Vector{Tv}, A::SparseMatrixCSC{Tv,Int64}, x::Vector{Tv}) where {Tv<:AbstractFloat}

    # check the size of matrix and vector
    if length(x) != A.n || length(y) != A.m
       error("dimension dismatch")
    end

    # Float32
    if eltype(y) == Float32
       ccall((:a_mul_b_ss_, spmatveclib),
             Void,
             (Ref{Int64}, Ref{Int64}, Ref{Float32}, Ref{Int64}, Ref{Int64}, Ref{Float32}, Ref{Float32}),
              A.m       , A.n       , A.nzval     , A.rowval  , A.colptr  , x           , y             )

    # Float64
    elseif eltype(y) == Float64
       ccall((:a_mul_b_rr_, spmatveclib),
             Void,
             (Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}),
              A.m       , A.n       , A.nzval     , A.rowval  , A.colptr  , x           , y             )
    end

    return nothing
end

"""
   y = A' * x, A is a sparse matrix in CSC format, x is a vector, y is updated in-place
can be used for both Float32 or Float64.
"""
function At_mul_b!(y::Vector{Tv}, A::SparseMatrixCSC{Tv,Ti}, x::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    # check the dimension
    if length(y) != A.n || length(x) != A.m
       error("dimension dismatch")
    end

    # Float32
    if eltype(x) == Float32
       ccall((:ac_mul_b_ss_, spmatveclib), Void,
             (Ref{Int64}, Ref{Int64}, Ref{Float32}, Ref{Int64}, Ref{Int64}, Ref{Float32}, Ref{Float32}),
              A.m       , A.n       , A.nzval     , A.rowval  , A.colptr  , x           , y             )

    # Float64
    elseif eltype(x) == Float64
       ccall((:ac_mul_b_rr_, spmatveclib), Void,
             (Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}),
              A.m       , A.n       , A.nzval     , A.rowval  , A.colptr  , x           , y             )
    end

    return nothing
end
