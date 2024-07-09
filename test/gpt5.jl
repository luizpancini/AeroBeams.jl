using SparseArrays
using LinearAlgebra

# Generate a square sparse matrix A using sprandn()
n = 5  # Size of the matrix
density = 0.5  # Density of the non-zero elements
A = sprandn(n, n, density)

# Ensure A is well-conditioned by adding a diagonal dominance
A += n * I

# Generate a sparse matrix B
B = sprandn(n, 1, density)

# Perform the factorization
F = lu(A)

# Solve the system using the factorization
X = F \ B

# Display the result
println("The solution X is:")
println(X)
