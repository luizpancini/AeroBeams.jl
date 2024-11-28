using AeroBeams, LinearAlgebra

# Beam rotation angle about a3-axis
α = π*30/180

# Beam rotation angle about a2'-axis
β = -π*30/180

# Hinge angle [rad]
θ = -π*45/180

# Wiener-Milenkovic parameters due to hinge angle, resolved in the local (b) basis
p_b = ypr_to_WM([0;θ;0])

# Wiener-Milenkovic parameters due to hinge angle, resolved in the global (A) basis (p_A = R0 * p_b)
p_A = rotation_tensor_E321([α;β;0]) * p_b

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],rotationParametrization="E321",p0=[α;β;0])

# Rotation constraints
rotationConstraints = Vector{RotationConstraint}()
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=hingeNode-1,slaveElementLocalID=hingeNode,masterDOF=1,slaveDOF=1,value=p_A[1],loadBalanceLocalNode=hingeNode+1))
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=hingeNode-1,slaveElementLocalID=hingeNode,masterDOF=2,slaveDOF=2,value=p_A[2],loadBalanceLocalNode=hingeNode+1))
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=hingeNode-1,slaveElementLocalID=hingeNode,masterDOF=3,slaveDOF=3,value=p_A[3],loadBalanceLocalNode=hingeNode+1))

# BCs
Fₕ = -1
Fₜ = -1
qᵢ = -1
qₒ = -1
q = (x1,t) -> ifelse.(x1.<=L*(hingeNode-1)/nElem,qᵢ,0) .+ ifelse.(x1.>L*(hingeNode-1)/nElem,qₒ,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
hingeForce = create_BC(name="hingeForce",beam=beam,node=hingeNode,types=["F3A"],values=[Fₕ])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[Fₜ])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
clampedHingedBeamRotated = create_Model(name="clampedHingedBeamRotated",beams=[beam],BCs=[clamp,hingeForce,tipForce],rotationConstraints=rotationConstraints)

# System solver
σ0 = 1
σstep = 0.5
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve the problem
problem = create_SteadyProblem(model=clampedHingedBeamRotated,systemSolver=NR)
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
p2_b = vcat([vcat(scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n1_b)[2],scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n2_b)[2]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Compute rotation, displacement and balance moment at the hinge node
Δp2_bHingeNode = 4*atand((p2_b[2*hingeNode+1]-p2_b[2*hingeNode-2])/4)
u3HingeNode = u3[2*hingeNode-1]
u3Tip = u3[end]
hingeBalanceM2 = -problem.model.rotationConstraints[1].balanceMoment
hingeBalanceM1 = -problem.model.rotationConstraints[2].balanceMoment
hingeBalanceM = sqrt(hingeBalanceM1^2+hingeBalanceM2^2)

# Display results
println("Rotation at hinge node = $Δp2_bHingeNode deg")
println("Displacement at hinge node = $u3HingeNode")
println("Displacement at tip = $u3Tip")
println("Moment necessary to impose the constraint at hinge node = $hingeBalanceM")

println("Finished clampedHingedBeamRotated.jl")