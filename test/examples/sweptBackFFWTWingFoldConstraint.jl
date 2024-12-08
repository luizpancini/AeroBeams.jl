using AeroBeams, LinearAlgebra

# Sweep angle [rad]
Λ = -30*π/180

# Flare angle [rad]
Γ = 30*π/180

# Fold angle [rad]
θ = -90*π/180

# Wiener-Milenkovic parameters due to hinge angle, resolved in the local (b) basis
p_b = ypr_to_WM([0; θ; 0])

# Wiener-Milenkovic parameters due to hinge angle, resolved in the global (A) basis (p_A = R0 * p_b)
R0 = rotation_tensor_E321([-Γ; 0; 0]) * rotation_tensor_E321([Λ; 0; 0])
p_A = R0 * p_b

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.5
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
nElem = 10
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix()],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],rotationParametrization="E321",p0=[Λ;0;0],aeroSurface=surf)

# Rotation constraints
rotationConstraints = Vector{RotationConstraint}()
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=hingeNode-1,slaveElementLocalID=hingeNode,masterDOF=1,slaveDOF=1,value=p_A[1],loadBalanceLocalNode=hingeNode+1))
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=hingeNode-1,slaveElementLocalID=hingeNode,masterDOF=2,slaveDOF=2,value=p_A[2],loadBalanceLocalNode=hingeNode+1))
push!(rotationConstraints,create_RotationConstraint(beam=beam,masterElementLocalID=hingeNode-1,slaveElementLocalID=hingeNode,masterDOF=3,slaveDOF=3,value=p_A[3],loadBalanceLocalNode=hingeNode+1))

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
sweptBackFFWTWingFoldConstraint = create_Model(name="sweptBackFFWTWingFoldConstraint",beams=[beam],BCs=[clamp],rotationConstraints=rotationConstraints)

# System solver
σ0 = 1
maxIter = 50
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve the problem
problem = create_SteadyProblem(model=sweptBackFFWTWingFoldConstraint,systemSolver=NR)
solve!(problem)

# Get nodal arclength positions, displacements and rotations over the beam
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[2],problem.nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
p3 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[3],problem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)

# Compute rotation at the hinge node
pOutboard = [p1[2*hingeNode-1]; p2[2*hingeNode-1]; p3[2*hingeNode-1]]
pInboard = [p1[2*hingeNode-2]; p2[2*hingeNode-2]; p3[2*hingeNode-2]]
ΔpHingeNode = pOutboard - pInboard
ypr = WM_to_ypr(ΔpHingeNode)*180/π

# Display results
println("Rotation at hinge node: yaw = $(ypr[1]) deg, pitch = $(ypr[2]) deg, roll = $(ypr[3])")

println("Finished sweptBackFFWTWingFoldConstraint.jl")