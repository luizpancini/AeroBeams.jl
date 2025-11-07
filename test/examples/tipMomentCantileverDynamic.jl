using AeroBeams, LinearAlgebra

# Beam
L = 1
EIy = 1
stiffnessMatrix = isotropic_stiffness_matrix(EIy=EIy)
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix])

# Time variables
tf = 1
Δt = tf/1e2

# BCs
M₀ = 2π*EIy/L
M = t -> M₀*t/tf
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipMoment = create_BC(name="tipMoment",beam=beam,node=nElem+1,types=["M2A"],values=[t->M(t)])

# Model
tipMomentCantileverDynamic = create_Model(name="tipMomentCantileverDynamic",beams=[beam],BCs=[clamp,tipMoment],units=create_UnitsSystem(length="in",force="lbf"))

# Create and solve the problem
problem = create_DynamicProblem(model=tipMomentCantileverDynamic,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
t = problem.savedTimeVector
κ2 = [problem.compElementalStatesOverTime[i][1].κ[2] for i in eachindex(t)]

println("Finished tipMomentCantileverDynamic.jl")