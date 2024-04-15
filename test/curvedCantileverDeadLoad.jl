using AeroBeams, LinearAlgebra, Plots

# Beam
R,θ = 100,π/4
A = 1
Iy = Iz = 1/12
J = Iy+Iz
L = R*θ
E,ν = 1e7,0
G = E/(2*(1+ν))
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*Iy,E*Iz])
nElem = 16
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],k=[0;0;1/R])

# BCs
F = 600
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["F3A"],values=[F])

# Model
curvedCantilever = Model(name="curvedCantilever",beams=[beam],BCs=[clamp,tipForce])

# Set system solver options
σ0 = 0.05
σstep = 0.05
NR = NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = SteadyProblem(model=curvedCantilever,systemSolver=NR)
solve!(problem)

# Get tip displacements  
tip_u1 = problem.nodalStatesOverσ[end][nElem].u_n2[1] 
tip_u2 = problem.nodalStatesOverσ[end][nElem].u_n2[2]
tip_u3 = problem.nodalStatesOverσ[end][nElem].u_n2[3] 

# Print 
println("Tip displacements:\nu1 = $tip_u1, u2 = $tip_u2, u3 = $tip_u3")

println("Finished curvedCantilever.jl")
