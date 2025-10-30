using AeroBeams, LinearAlgebra

# Beam
L = 1
EIy = 1
stiffnessMatrix = isotropic_stiffness_matrix(EIy=EIy)
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipMoment = create_BC(name="tipMoment",beam=beam,node=nElem+1,types=["M2A"],values=[2*π*EIy/L])

# Model
tipMomentCantilever = create_Model(name="tipMomentCantilever",beams=[beam],BCs=[clamp,tipMoment],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tipMomentCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][nElem].θ_n2 for i in 1:length(σVector)]

println("Finished tipMomentCantilever.jl")