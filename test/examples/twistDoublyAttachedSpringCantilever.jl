using AeroBeams, LinearAlgebra

# Beam 
L = 1
GJ = 1
nElem = 10
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(GJ=GJ)])

# Spring (attached at middle and tip nodes)
k = 1
spring = create_Spring(elementsIDs=[midElem,nElem],nodesSides=[2,2],kp=[k,0,0])
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
M₀ = 1e-1
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
moment = create_BC(name="moment",beam=beam,node=nElem+1,types=["M1A"],values=[M₀])

# Model
twistDoublyAttachedSpringCantilever = create_Model(name="twistDoublyAttachedSpringCantilever",beams=[beam],BCs=[clamp,moment])

# System solver
σ0 = 1
maxIter = 20
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=twistDoublyAttachedSpringCantilever,systemSolver=NR)
solve!(problem)

# Get outputs
x1_n = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
M1 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[1],problem.nodalStatesOverσ[end][e].M_n2[1]) for e in 1:nElem]...)

# Analytical solution
pTipAnalytical = M₀*L * (1/(2GJ) + 1/(2GJ+k*L))
MsAnalytical = M₀*(1-1/(1+k*L/(2GJ)))

# Relative errors
ϵp = p1[end]/pTipAnalytical - 1
ϵMs = norm(spring.Ms)/MsAnalytical - 1

# Display results
println("Relative error in tip twist rotation = $ϵp")
println("Relative error in spring moment = $ϵMs")

println("Finished twistDoublyAttachedSpringCantilever.jl")