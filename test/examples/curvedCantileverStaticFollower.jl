using AeroBeams, LinearAlgebra, DelimitedFiles

# Beam
R,θ = 100,π/4
L = R*θ
A,Iy,Iz,J = 1,1/12,1/12,1/6
E,G = 1e7,5e6
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*Iy,E*Iz])
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],k=[0;0;1/R])

# BCs
F = 3000
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[F])

# Model
curvedCantileverStaticFollower = create_Model(name="curvedCantileverStaticFollower",beams=[beam],BCs=[clamp,tipForce])

# Set system solver options
σ0 = 0.0
σstep = 1e-2
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=curvedCantileverStaticFollower,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u2 = [problem.nodalStatesOverσ[i][nElem].u_n2[2] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]

# Reference solution by Simo and Vu-Quoc (1986)
u1_ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "curvedCantileverStaticFollower", "u1.txt"))
u2_ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "curvedCantileverStaticFollower", "u2.txt"))
u3_ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "curvedCantileverStaticFollower", "u3.txt"))

println("Finished curvedCantileverStaticFollower.jl")
