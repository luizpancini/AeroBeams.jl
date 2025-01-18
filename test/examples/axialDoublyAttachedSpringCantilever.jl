using AeroBeams, LinearAlgebra

# Beam 
L = 1
EA = 1
nElem = 10
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA)])

# Spring (attached at middle and tip nodes)
k = 2
spring = create_Spring(elementsIDs=[midElem,nElem],nodesSides=[2,2],ku=[k;0;0])
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
F₀ = 1
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
force = create_BC(name="force",beam=beam,node=nElem+1,types=["F1A"],values=[F₀])

# Model
axialDoublyAttachedSpringCantilever = create_Model(name="axialDoublyAttachedSpringCantilever",beams=[beam],BCs=[clamp,force])

# System solver
σ0 = 1
maxIter = 20
relTol = 1e-9
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=axialDoublyAttachedSpringCantilever,systemSolver=NR)
solve!(problem)

# Get outputs
x1_n = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
F1 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[1],problem.nodalStatesOverσ[end][e].F_n2[1]) for e in 1:nElem]...)

# Analytical solution
uTipAnalytical = F₀*L * (1/(2EA) + 1/(2EA+k*L))
FsAnalytical = F₀*(1-1/(1+k*L/(2EA)))

# Relative errors
ϵu = u1[end]/uTipAnalytical - 1
ϵFs = norm(spring.Fs)/FsAnalytical - 1

# Display results
println("Relative error in tip axial displacement = $ϵu")
println("Relative error in spring force = $ϵFs")

println("Finished axialDoublyAttachedSpringCantilever.jl")