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
u₀ = 1/3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipDisp = create_BC(name="tipDisp",beam=beam,node=nElem+1,types=["u1A"],values=[u₀])

# Model
axialDoublyAttachedSpringCantilever2 = create_Model(name="axialDoublyAttachedSpringCantilever2",beams=[beam],BCs=[clamp,tipDisp])

# System solver
σ0 = 1
maxIter = 20
relTol = 1e-9
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=axialDoublyAttachedSpringCantilever2,systemSolver=NR)
solve!(problem)

# Get outputs
x1_n = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
F1 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[1],problem.nodalStatesOverσ[end][e].F_n2[1]) for e in 1:nElem]...)

# Analytical solution
FrootAnalytical = u₀/L / (1/(2EA) + 1/(2EA+k*L))
FsAnalytical = k * (FrootAnalytical/(1+k*L/(2EA)) * L/(2EA))

# Relative errors
ϵFr = F1[1]/FrootAnalytical - 1
ϵFs = norm(spring.Fs)/FsAnalytical - 1

# Display results
println("Relative error in root reaction = $ϵFr")
println("Relative error in spring force = $ϵFs")

println("Finished axialDoublyAttachedSpringCantilever2.jl")