using AeroBeams, LinearAlgebra

# Angles of rotation in the Euler 3-2-1 parametrization:
# Hinge angle about a3 (IP)
α = 0

# Hinge angle about a2 (OOP)
β = -π*45/180

# Hinge angle about a1 (twist)
γ = -π*20/180

# Equivalent Wiener-Milenkovic rotation parameters
p_WM = ypr_to_WM([α,β,γ])

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.25
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
EIy = 1
GJ = 1e1
EIz = 1e2
nElem = 20
hingeNode = 15
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(GJ=GJ,EIy=EIy,EIz=EIz)],hingedNodes=[hingeNode],hingedNodesDoF=[trues(3)],aeroSurface=surf)

# Rotation constraints
rotationConstraints = Vector{RotationConstraint}()
inboardElem = hingeNode-1
outboardElem = hingeNode
loadBalanceLocalNode = hingeNode+1
# γ
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,masterDOF=1,slaveDOF=1,value=p_WM[1],loadBalanceLocalNode=loadBalanceLocalNode))
# β
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,masterDOF=2,slaveDOF=2,value=p_WM[2],loadBalanceLocalNode=loadBalanceLocalNode))
# α
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,masterDOF=3,slaveDOF=3,value=p_WM[3],loadBalanceLocalNode=loadBalanceLocalNode))

# BCs
Fₕ = -0
Fₜ = -0.5
qᵢ = 2
qₒ = 2
q = (x1,t) -> ifelse.(x1.<=L*(hingeNode-1)/nElem,qᵢ,0) .+ ifelse.(x1.>L*(hingeNode-1)/nElem,qₒ,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; -q(x1,t); q(x1,t)]])
hingeForce = create_BC(name="hingeForce",beam=beam,node=hingeNode,types=["F3A"],values=[Fₕ])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[Fₜ])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
clampedUniversalHingeBeam = create_Model(name="clampedUniversalHingeBeam",beams=[beam],BCs=[clamp,hingeForce,tipForce],rotationConstraints=rotationConstraints)

# System solver
σ0 = 1
σstep = 0.5
maxIter = 100
relTol = 1e-9
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve the problem
problem = create_SteadyProblem(model=clampedUniversalHingeBeam,systemSolver=NR)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[2],problem.nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p1 = vcat([vcat(scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n1)[1],scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n2)[1]) for e in 1:nElem]...)
p2 = vcat([vcat(scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n1)[2],scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n2)[2]) for e in 1:nElem]...)
p3 = vcat([vcat(scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n1)[3],scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n2)[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Compute rotations, displacement and balance moment at the hinge node
Δp_WMHingeNode = [p1[2*hingeNode]-p1[2*hingeNode-2]; p2[2*hingeNode]-p2[2*hingeNode-2]; p3[2*hingeNode]-p3[2*hingeNode-2]]
R,_ = rotation_tensor_WM(Δp_WMHingeNode)
yaw,pitch,roll = ypr_from_rotation_tensor(R)
ΔpHingeNode = [yaw,pitch,roll]*180/π
u3HingeNode = u3[2*hingeNode-1]
hingeBalanceM1 = -problem.model.rotationConstraints[1].balanceMoment
hingeBalanceM2 = -problem.model.rotationConstraints[2].balanceMoment
hingeBalanceM3 = -problem.model.rotationConstraints[3].balanceMoment
hingeBalanceM = sqrt(hingeBalanceM1^2+hingeBalanceM2^2+hingeBalanceM3^2)

# Display results
println("Twist rotation at hinge node = $(ΔpHingeNode[3]) deg")
println("OOP rotation at hinge node = $(ΔpHingeNode[2]) deg")
println("IP rotation at hinge node = $(ΔpHingeNode[1]) deg")
println("Displacement at hinge node = $u3HingeNode")
println("Moment necessary to impose the hinge angle = $hingeBalanceM")

println("Finished clampedUniversalHingeBeam.jl")