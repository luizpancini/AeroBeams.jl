using AeroBeams, LinearAlgebra

# Hinge flare angle
Λ = 30*π/180

# Beam rotation angle about axis a3
α = 0*π/180

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.5
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Gravity
g = 9.80665

# Beam
L = 1
EIy = 1e1
ρA = 1e-1
nElem = 10
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],aeroSurface=surf,rotationParametrization="E321",p0=[α;0;0])

# Hinge axis constraint (with hinge axis defined in the local, undeformed beam basis)
pHValue = 4*tan(pi/2/4)
solutionMethod = "addedResidual"
updateAllDOFinResidual = false
localHingeAxis = rotation_tensor_E321([-Λ; 0; 0]) * [0; 1; 0]
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,localHingeAxis=localHingeAxis,pHValue=pHValue)

 # Spring around hinge
kT = kOOP = kIP = 1e-6
spring = create_Spring(elementsIDs=[hingeNode-1,hingeNode],nodesSides=[1,2],kp=[kT,kOOP,kIP])
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# System solver
σ0 = 1
maxIter = 50
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Model
clampedFlaredHingedBeam = create_Model(name="clampedFlaredHingedBeam",beams=[beam],BCs=[clamp],gravityVector=[0,0,-g],hingeAxisConstraints=[hingeAxisConstraint])

# Create and solve the problem
problem = create_SteadyProblem(model=clampedFlaredHingedBeam,systemSolver=NR)
solve!(problem)

# Get outputs
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[2],problem.nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
p3 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[3],problem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Get nodal arclength positions
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Extract rotation angle and hinge moment about the hinge node
ϕHinge = problem.model.hingeAxisConstraints[1].ϕ*180/π
hingeMoment = problem.model.hingeAxisConstraints[1].hingeMoment

# Display results
println("Rotation at hinge node: $ϕHinge deg")
println("Hinge moment: $hingeMoment N.m")

println("Finished clampedFlaredHingedBeam.jl")