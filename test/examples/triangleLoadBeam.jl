using AeroBeams, LinearAlgebra

# Beam
L = 1.0
EIy = 1e7
∞ = 1e12
stiffnessMatrix = diagm([∞,∞,∞,∞,EIy,∞])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
q₀ = 1e6
q = (x1,t) -> q₀*x1
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
roller = create_BC(name="roller",beam=beam,node=1,types=["u3A"],values=[0])

# Model
triangleLoadBeam = create_Model(name="triangleLoadBeam",beams=[beam],BCs=[clamp,roller])

# Create and solve the problem
problem = create_SteadyProblem(model=triangleLoadBeam,getLinearSolution=true)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

println("Finished triangleLoadBeam.jl")