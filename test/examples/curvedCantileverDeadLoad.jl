using AeroBeams, LinearAlgebra

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
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],k=[0;0;1/R])

# BCs
F = 600
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["F3A"],values=[F])

# Model
curvedCantilever = create_Model(name="curvedCantilever",beams=[beam],BCs=[clamp,tipForce])

# Set system solver options
σ0 = 0.05
σstep = 0.05
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=curvedCantilever,systemSolver=NR)
solve!(problem)

# Get tip displacements  
tip_u1 = problem.nodalStatesOverσ[end][nElem].u_n2[1] 
tip_u2 = problem.nodalStatesOverσ[end][nElem].u_n2[2]
tip_u3 = problem.nodalStatesOverσ[end][nElem].u_n2[3] 

# Reference solution by Bathe & Bolourchi - Large displacement analysis of three-dimensional beam structures (1979)
u1Ref = -23.5
u2Ref = -13.4
u3Ref = 53.4

# Compute and print relative errors
ϵu1 = 1 - tip_u1/u1Ref
ϵu2 = 1 - tip_u2/u2Ref
ϵu3 = 1 - tip_u3/u3Ref
println("Relative errors: u1 = $ϵu1, u2 = $ϵu2, u3 = $ϵu3")

println("Finished curvedCantileverDeadLoad.jl")
