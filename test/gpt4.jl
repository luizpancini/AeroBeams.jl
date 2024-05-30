using LinearAlgebra

# Define the nonlinear system of equations
function F(x)
    return [x[1]^2 + x[2]^2 - 4, x[1]*x[2] - 1]
end

# Define the Jacobian of the system
function J(x)
    return [2*x[1] 2*x[2]; x[2] x[1]]
end

# Regularized inverse function
function regularized_inverse(J, λ=1e-6)
    n = size(J, 1)
    return inv(J + λ * I(n))
end

# Line search function
function line_search(f, x, p, Jx, c1=1e-4, ρ=0.5)
    fx = norm(f(x))^2
    β = dot(p, Jx * p)
    α = 1.0
    while norm(f(x + α * p))^2 > fx + c1 * α * β
        α *= ρ
    end
    return α
end

function newton_solver(F, J, x0; tol=1e-6, max_iter=100, λ=1e-8)
    x = x0
    for k in 1:max_iter
        Fx = F(x)
        Jx = J(x)
        
        # Check for singularity and regularize if needed
        if abs(det(Jx)) < tol
            Jx_inv = regularized_inverse(Jx, λ)
        else
            Jx_inv = inv(Jx)
        end
        
        # Solve for the Newton step
        p = -Jx_inv * Fx
        
        # Perform line search
        α = line_search(F, x, p, Jx)
        
        # Update the solution
        x += α * p
        
        # Check for convergence
        if norm(F(x)) < tol
            println("Converged in $k iterations")
            return x
        end
    end
    println("Did not converge")
    return x
end

# Initial guess
x0 = [1.0, 1.0]

# Solve the system
# solution = newton_solver(F, J, x0)

solver = (x0) -> newton_solver(F, J, x0)
solution = solver(x0)

println("Solution: ", solution)
println("F(solution): ", F(solution))
