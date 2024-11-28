using AeroBeams, LinearAlgebra, DelimitedFiles

# Beam 
R,θ = 0.5,π
L = R*θ
A,Iy = 1e-4,0.5e-8
E = 72e9
∞ = 1e12
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-π/2;0],k=[0;1/R;0])

# BCs
F = -5.5e3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,types=["Ff1A"],values=[F])

# Model
transverselyForcedArch = create_Model(name="transverselyForcedArch",beams=[beam],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=transverselyForcedArch,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][end].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][end].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][end].θ_n2 for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "transverselyForcedArch", "u1.txt"))
u3Ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "transverselyForcedArch", "u3.txt"))
θRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "transverselyForcedArch", "theta.txt"))

println("Finished transverselyForcedArch.jl")