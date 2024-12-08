using AeroBeams, LinearAlgebra

# Hinge flare angle
Λ = 30*π/180

# Beam rotation angle about axis a3
α = -30*π/180

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.5
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
EIy = 1
nElem = 16
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],aeroSurface=surf,rotationParametrization="E321",p0=[α;0;0])

# Hinge axis constraint (defined in the local, undeformed beam basis)
localHingeAxis = rotation_tensor_E321([-Λ; 0; 0]) * [0; 1; 0]
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,elementLocalID=hingeNode,localHingeAxis=localHingeAxis,loadBalanceLocalNode=hingeNode+1,foldGuessValue=π/2)

# Set spring around hinge (needed to avoid singular Jacobian)
k = 1e-4
doublyAttachedSpring = create_Spring(basis="A",elementsIDs=[hingeNode-1,hingeNode],nodesSides=[1,2],kTwist=k*1e0,kIPBending=k*1e0,kOOPBending=k*1e0)
add_spring_to_beams!(beams=[beam,beam],spring=doublyAttachedSpring)

# BCs
Fᵢ = 1e-2
Fₜ = 1e-2
q₀ = -1e-1
q = (x1,t) -> q₀
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
tipInplaneForce = create_BC(name="tipInplaneForce",beam=beam,node=nElem+1,types=["Ff2A"],values=[-Fᵢ])
tipNormalForce = create_BC(name="tipNormalForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[Fₜ])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# System solver
σ0 = 1
maxIter = 50
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Model
clampedHingedSpringedBeam = create_Model(name="clampedHingedSpringedBeam",beams=[beam],BCs=[clamp,tipNormalForce,tipInplaneForce],hingeAxisConstraints=[hingeAxisConstraint])

# Create and solve the problem
problem = create_SteadyProblem(model=clampedHingedSpringedBeam,systemSolver=NR)
solve!(problem)

# Get outputs
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[2],problem.nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p1_b = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1_b[1],problem.nodalStatesOverσ[end][e].p_n2_b[1]) for e in 1:nElem]...)
p2_b = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1_b[2],problem.nodalStatesOverσ[end][e].p_n2_b[2]) for e in 1:nElem]...)
p3_b = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1_b[3],problem.nodalStatesOverσ[end][e].p_n2_b[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
ps2_b = vcat([vcat(scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n1_b)[2],scaled_rotation_parameters(problem.nodalStatesOverσ[end][e].p_n2_b)[2]) for e in 1:nElem]...)

# Get nodal arclength positions
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Compute Euler angles at the hinge node
pOutboard = [p1_b[2*hingeNode-1]; p2_b[2*hingeNode-1]; p3_b[2*hingeNode-1]]
pInboard = [p1_b[2*hingeNode-2]; p2_b[2*hingeNode-2]; p3_b[2*hingeNode-2]]
ΔpHingeNode = pOutboard - pInboard
ypr = WM_to_ypr(ΔpHingeNode)*180/π

# Display results
println("Rotation at hinge node: yaw = $(ypr[1]) deg, pitch = $(ypr[2]) deg, roll = $(ypr[3])")

println("Finished clampedHingedSpringedBeam.jl")