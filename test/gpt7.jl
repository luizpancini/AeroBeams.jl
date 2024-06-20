using LinearAlgebra

# Example matrices
A = [1 2 3; 4 5 6; 7 8 9]  # 3x3 square matrix
B = [1 2 3; 4 5 6]          # 2x3 non-square matrix

# Check if the matrices are square
println(issquare(A))  # Output: true
println(issquare(B))  # Output: false
