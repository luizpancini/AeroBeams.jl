using AeroBeams, LinearAlgebra

# Beam 1
L1 = 2
EIy1 = 1e16
nElem1 = 10
beam1 = create_Beam(name="beam1",length=L1,nElements=nElem1,S=[isotropic_stiffness_matrix(EIy=EIy1)])

# Beam 2
L2 = 1
EIy2 = 1e4
nElem2 = 20
beam2 = create_Beam(name="beam2",length=L2,nElements=nElem2,S=[isotropic_stiffness_matrix(EIy=EIy2)],rotationParametrization="E321",p0=[0;π/2;0])

# Spring
kSpring = 1e4
spring = create_Spring(basis="b",elementsIDs=[div(nElem1,2),nElem2],nodesSides=[2,2],ku=[kSpring,kSpring,kSpring])
add_spring_to_beams!(beams=[beam1,beam2],spring=spring)

# BCs
F = 1e2
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam2,node=nElem2+1,types=["F3b"],values=[F])

# Model
springedLFrame = create_Model(name="springedLFrame",beams=[beam1,beam2],BCs=[clamp,tipForce])

# System solver
maxIter = 50
NR = create_NewtonRaphson(maximumIterations=maxIter,displayStatus=true)

# Create and solve the problem
problem = create_SteadyProblem(model=springedLFrame,systemSolver=NR)
solve!(problem)

# Get nodal arclength positions and displacement over tip beam
x1 = vcat([vcat(problem.model.beams[2].elements[e].x1_n1,problem.model.beams[2].elements[e].x1_n2) for e in 1:nElem2]...)
u3_b = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1_b[3],problem.nodalStatesOverσ[end][e].u_n2_b[3]) for e in nElem1+1:nElem1+nElem2]...)

# Tip u3 displacement (in local basis)
println("Tip disp = $(u3_b[end])")

println("Finished springedLFrame.jl")