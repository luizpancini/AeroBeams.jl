using SparseArrays, LinearAlgebra

# Define sparse matrices A and B
A::SparseMatrixCSC{Float64,Int64} = sparse([1, 2, 3], [1, 2, 3], [4.0, 5.0, 6.0])
B::SparseMatrixCSC{Float64,Int64} = sparse([1, 2, 3], [1, 2, 3], [7.0, 8.0, 9.0])
A::SparseMatrixCSC{Float64,Int64} = sprandn(3,3,0.99)
B::SparseMatrixCSC{Float64,Int64} = sprandn(3,3,0.5)

# Compute A \ B
X = A \ Array(B)

# Display the result
println("The solution X is:")
println(X)
