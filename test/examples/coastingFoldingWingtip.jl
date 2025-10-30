using AeroBeams

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.25
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeElemNorm = 3/4
hingeElem = floor(Int,hingeElemNorm*nElem)
hingedNode = hingeElem+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingedNode],hingedNodesDoF=[[true,true,true]],aeroSurface=surf)

# Hinge flare angle and local axis
flareAngle = 20*π/180
localHingeAxis = rotation_tensor_E321([-flareAngle; 0; 0]) * AeroBeams.a2

# Solution method for hinge constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = false

# Hinge axis
hingeAxisConstraint = create_HingeAxisConstraint(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,beam=beam,localHingeAxis=localHingeAxis)

# Spring
kSpring = 1e-4
spring = create_Spring(elementsIDs=[hingeElem,hingeElem+1],nodesSides=[1,2],kp=[0,kSpring,0])
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
qᵢ = -1e1
q₀ = -1e-1
q = (x1,t) -> ifelse.(x1.<L*hingeElemNorm,qᵢ,0) .+ ifelse.(x1.>=L*hingeElemNorm,q₀,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
coastingFoldingWingtip = create_Model(name="coastingFoldingWingtip",beams=[beam],BCs=[clamp],hingeAxisConstraints=[hingeAxisConstraint])

# System solver
σ0 = 1
maxIter = 100
desIter = 50
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,desiredIterations=desIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=coastingFoldingWingtip,systemSolver=NR)
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

# Display results
println("Coast angle = $ϕHinge deg")

println("Finished coastingFoldingWingtip.jl")