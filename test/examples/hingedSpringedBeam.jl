using AeroBeams, LinearAlgebra

# Option to solve linear problem
linear = false

# Beam 
L = 2
EA = 1e6
EIy = 1e6
nElem = 100
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(EA=EA,EIy=EIy)],hingedNodes=[midElem+1],hingedNodesDoF=[[false,true,false]])

# Spring
kOOPBending = 1e5
spring = create_Spring(basis="A",elementsIDs=[midElem,midElem+1],nodesSides=[1,2],kOOPBending=kOOPBending)
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
q₀ = 1e4
q = (x1,t) -> ifelse.(x1.>=L/2,-q₀,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
hingedSpringedBeam = create_Model(name="hingedSpringedBeam",beams=[beam],BCs=[pin,clamp])

# System solver
maxIter = 1_000
relTol = 1e-8
absTol = 1e-8
NR = create_NewtonRaphson(maximumIterations=maxIter,absoluteTolerance=absTol,relativeTolerance=relTol,displayStatus=true)

# Create and solve the problem
problem = create_SteadyProblem(model=hingedSpringedBeam,systemSolver=NR,getLinearSolution=linear)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Display midpoint displacement
println("Midpoint u3 = $(u3[midElem+1])")

println("Finished hingedSpringedBeam.jl")