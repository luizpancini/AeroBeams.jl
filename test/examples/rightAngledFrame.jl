using AeroBeams, LinearAlgebra, DelimitedFiles

# Beam frame
L = 0.24
A,Iy = 0.18e-4,0.135e-8
E = 7.124e10
ν = 0.3
G = E/(2*(1+ν))
∞ = 1e12
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[stiffnessMatrix])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[stiffnessMatrix],rotationParametrization="E321",p0=[0;π/2;0])

# BCs
F = 5e3
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam2,node=nElem+1,types=["Ff3b"],values=[-F])

# Model
rightAngledFrame = create_Model(name="rightAngledFrame",beams=[beam1,beam2],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=rightAngledFrame,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][end].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][end].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][end].θ_n2 for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "rightAngledFrame", "u1.txt"))
u3Ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "rightAngledFrame", "u3.txt"))
θRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "rightAngledFrame", "theta.txt"))

println("Finished rightAngledFrame.jl")