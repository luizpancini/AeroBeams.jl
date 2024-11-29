using AeroBeams, LinearAlgebra

# Option for linear solution
linear = true

# Beam
a = b = 1 
L = a+b
EIy = 1
nElem = 20
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[midElem+1],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = -1
F₀ = -1
q = (x1,t) -> ifelse.(x1.>=L/2,q₀,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
force = create_BC(name="force",beam=beam,node=midElem+1,types=["F3A"],values=[F₀])
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
biclampedHingedBeam = create_Model(name="biclampedHingedBeam",beams=[beam],BCs=[clamp1,clamp2,force])

# System solver (for nonlinear solution)
NR = create_NewtonRaphson(displayStatus=false)

# Create and solve the problem
problem = create_SteadyProblem(model=biclampedHingedBeam,getLinearSolution=linear,systemSolver=NR)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Linear FEM solution (problem 5.20 of REDDY - An Introduction to the Finite Element Method - 3rd Ed. - 2006)
u3MidRef = 1/EIy*a^3/(a^3+b^3)*(q₀*b^4/8+F₀*b^3/3)
p2RightRef = 1/EIy*b^3/(a^3+b^3)*(q₀*(8*a^3-b^3)/48+F₀*a^3/(2*b))

# The solution for the angle of rotation at the left of the hinge is not given by Reddy, but it is close to the following (verified using a linear beam FE code):
p2LeftRef = -1/EIy*b^3/(a^3+b^3)*(q₀*(8*a^3+b^3)/48+F₀*a^3/(2*b))

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

println("Finished biclampedHingedBeam.jl")