using AeroBeams

# Beam
L = 1
EIy = 1e4
θ = 30*π/180
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],rotationParametrization="E321",p0=[0;θ;0])

# Spring
ku = 1e4*[0; 0; 1]
spring = create_Spring(basis="b",elementsIDs=[nElem],nodesSides=[2],ku=ku)
add_springs_to_beam!(beam=beam,springs=[spring])

# BCs
F = 1e2
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["F3b"],values=[F])

# Model
cantileverWithTipSpring = create_Model(name="cantileverWithTipSpring",beams=[beam],BCs=[clamp,tipForce])

# Create and solve the problem
problem = create_SteadyProblem(model=cantileverWithTipSpring)
solve!(problem)

# Get nodal arclength positions, displacement and loads over the beam
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1_b[3],problem.nodalStatesOverσ[end][e].u_n2_b[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Tip displacement
println("Tip disp = $(u3[end])")

println("Finished cantileverWithTipSpring.jl")