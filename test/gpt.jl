using ForwardDiff, FiniteDifferences

function f(x)
    p,u = x[1:3],x[4]
    # Example function
    y1 = p[1] * p[2] + sin(p[3])
    y2 = exp(p[1]) - p[2]^u
    return [y1, y2]
end

# Define inputs
p = [1.0, 2.0, 3.0]
u = 2
x = vcat([p,u]...)

# Calculate the Jacobian matrix of function f at input vector x
jacobian_matrix_AD = ForwardDiff.jacobian(f, x)

# Print the Jacobian matrix
println("Jacobian matrix using ForwardDiff at x = $x:")
println(jacobian_matrix_AD)

# Calculate the Jacobian matrix of function f at input vector x
jacobian_matrix_FD = FiniteDifferences.jacobian(central_fdm(5, 1), f, x)

# Print the Jacobian matrix
println("Jacobian matrix using FiniteDifferences at x = $x:")
println(jacobian_matrix_FD)
