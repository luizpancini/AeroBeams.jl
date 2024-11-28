using AeroBeams, LinearAlgebra

# Option for linear solution
linear = true

# Beam 
L = 2
EIy = 1
nElem = 20
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[midElem+1],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = -1
q = (x1,t) -> ifelse.(x1.>=L/2,q₀,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
roller = create_BC(name="roller",beam=beam,node=1,types=["u3A"],values=[0])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
hingedBeam = create_Model(name="hingedBeam",beams=[beam],BCs=[roller,clamp])

# Create and solve the problem
problem = create_SteadyProblem(model=hingedBeam,getLinearSolution=linear)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beam
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Linear FEM solution (example 5.2.4 of REDDY - An Introduction to the Finite Element Method - 3rd Ed. - 2006)
u3MidRef = q₀*(L/2)^4/(8*EIy)
p2LeftRef = -q₀*(L/2)^4/(8*L/2*EIy)
p2RightRef = q₀*(L/2)^3/(6*EIy)

# AeroBeams solution
u3Mid = u3[2*midElem+1]
p2Left = p2[2*midElem]
p2Right = p2[2*midElem+1]

# Compute relative errors
ϵ_u3Mid = 1 - u3Mid/u3MidRef
ϵ_p2Left = 1 - p2Left/p2LeftRef
ϵ_p2Right = 1 - p2Right/p2RightRef

# Display relative errors
println("Relative errors: u3Mid = $ϵ_u3Mid, p2Left = $ϵ_p2Left, p2Right = $ϵ_p2Right")

println("Finished hingedBeam.jl")