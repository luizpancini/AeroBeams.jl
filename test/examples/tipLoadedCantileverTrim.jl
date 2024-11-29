using AeroBeams, LinearAlgebra

# Beam
L = 1
EI = 333.333
∞ = 1e14
stiffnessMatrix = diagm([∞,∞,∞,∞,EI,∞])
nElem = 30
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix])

# BCs
F = 1
trimF3Guess,trimM2Guess = 0,0
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["F3A"],values=[-F])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clampReactions = create_BC(name="clampReactions",beam=beam,node=1,types=["F3A","M2A"],values=[trimF3Guess,trimM2Guess],toBeTrimmed=[true,true])

# Model
tipLoadedCantileverTrim = create_Model(name="tipLoadedCantileverTrim",beams=[beam],BCs=[clamp,tipForce,clampReactions])

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=50,displayStatus=false)

# Create and solve the problem
problem = create_TrimProblem(model=tipLoadedCantileverTrim,systemSolver=NR)
solve!(problem)

# Get numerical solution 
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Analytical solution (valid for small displacements)
u3_analytical = x1 -> -F/(6*EI)*x1.^2 .* (3*L .- x1)
F3_analytical = x1 -> -F
M2_analytical = x1 -> F*L * (1 .- x1)

println("Finished tipLoadedCantileverTrim.jl")