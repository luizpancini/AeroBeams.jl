using AeroBeams, LinearAlgebra, DelimitedFiles

# Beam
L = 100
E = 420e6
A,Iy = 1,1/12
∞ = 1e14
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*Iy,∞])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
F = 130e3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[F])

# Model
tipFollowerForceCantilever = create_Model(name="tipFollowerForceCantilever",beams=[beam],BCs=[clamp,tipFollowerForce],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tipFollowerForceCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "tipFollowerForceCantilever", "u1.txt"))
u3Ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "tipFollowerForceCantilever", "u3.txt"))

println("Finished tipFollowerForceCantilever.jl")