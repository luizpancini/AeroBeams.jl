using AeroBeams

# Beam 
L = 1
EIy = 1
nElem = 20
midElem = div(nElem,2)
hingedNode = midElem+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingedNode],hingedNodesDoF=[[true,true,true]])

# Hinge angle
hingeAngle = 90*π/180

# Solution method for hinge constraint
solutionMethod = "appliedMoment"

# Hinge constraint
hingeAxisConstraint = create_HingeAxisConstraint(solutionMethod=solutionMethod,beam=beam,localHingeAxis=AeroBeams.a2,pHValue=4*tan(hingeAngle/4))

# Spring
kSpring = hingeAngle/(4*tan(hingeAngle/4))
spring = create_Spring(elementsIDs=[midElem,midElem+1],nodesSides=[1,2],kp=[0,kSpring,0])
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
drivenSpringedHingedBeam = create_Model(name="drivenSpringedHingedBeam",beams=[beam],BCs=[clamp],hingeAxisConstraints=[hingeAxisConstraint])

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-9
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=drivenSpringedHingedBeam,systemSolver=NR)
solve!(problem)

# Get outputs
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
x1_e = [problem.model.beams[1].elements[e].x1 for e in 1:nElem]
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[2],problem.nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
p3 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[3],problem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Get rotation at the hinge node
pHinge = problem.model.hingeAxisConstraints[1].pH
ϕHinge = problem.model.hingeAxisConstraints[1].ϕ*180/π

# Analytical root moment and relative error
M2rootAnalytical = hingeAngle
ϵM2root = 1 - M2[1]/M2rootAnalytical

# Check results
println("Hinge angle = $ϕHinge deg")
println("Relative error on root bending moment = $ϵM2root")

println("Finished drivenSpringedHingedBeam.jl")