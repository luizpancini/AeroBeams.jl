using AeroBeams, LinearAlgebra

# Method for imposition of constraint ("rotationConstraint" or "hingeAxisConstraint")
foldConstraintMethod = "rotationConstraint"

# Sweep angle [rad]
Λ = -30*π/180

# Flare angle [rad]
Γ = 30*π/180

# Fold angle [rad]
θ = -90*π/180

# Flag to update hinge axis for computation of rotation
updateHingeAxis = false

# Hinge axis (defined in the local, undeformed beam basis)
localHingeAxis = rotation_tensor_E321([-Γ; 0; 0]) * AeroBeams.a2

# Set the reference DOF as the greatest component of hinge axis
refDOF = argmax(abs.(localHingeAxis))

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.25
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = 3*div(nElem,4)+1
outboardElem = hingeNode
inboardElem = outboardElem-1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],rotationParametrization="E321",p0=[Λ;0;0],aeroSurface=surf)

# Hinge axis constraint
ΔpValue = foldConstraintMethod == "hingeAxisConstraint" ? 4*tan(θ/4) : nothing
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,localHingeAxis=localHingeAxis,loadBalanceLocalNode=hingeNode+1,ΔpValue=ΔpValue,updateHingeAxis=updateHingeAxis)

# Fold angle constraint (if applicable)
rotationConstraints = Vector{RotationConstraint}()
if foldConstraintMethod == "rotationConstraint"
    foldAngleConstraint = create_RotationConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,masterDOF=refDOF,slaveDOF=refDOF,value=4*tan(θ/4)*hingeAxisConstraint.initialHingeAxis[refDOF],loadBalanceLocalNode=hingeNode+1)
    push!(rotationConstraints,foldAngleConstraint)
end

# BCs
Fₜ = -1e-1
qᵢ = -20e-0
qₒ = -0e-1
q = (x1,t) -> ifelse.(x1.<=L*(hingeNode-1)/nElem,qᵢ,0) .+ ifelse.(x1.>L*(hingeNode-1)/nElem,qₒ,0)
add_loads_to_beam!(beam,loadTypes=["ff_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[Fₜ])

# Model
sweptBackFFWTWingHingeAxisAndFoldConstraint = create_Model(name="sweptBackFFWTWingHingeAxisAndFoldConstraint",beams=[beam],BCs=[clamp,tipFollowerForce],hingeAxisConstraints=[hingeAxisConstraint],rotationConstraints=rotationConstraints)

# System solver
σ0 = 0.5
maxIter = 100
relTol = 1e-8
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Create and solve the problem
problem = create_SteadyProblem(model=sweptBackFFWTWingHingeAxisAndFoldConstraint,systemSolver=NR)
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
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
p3 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[3],problem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Get rotation angle and balance moment at the hinge node
ΔϕHinge = problem.model.hingeAxisConstraints[1].Δϕ*180/π
hingeBalanceM = foldConstraintMethod == "hingeAxisConstraint" ? norm(problem.model.hingeAxisConstraints[1].balanceMoment) : norm(vcat(problem.model.rotationConstraints[1].balanceMoment,problem.model.hingeAxisConstraints[1].balanceMoment))

# Display results
println("Hinge rotation angle = $ΔϕHinge deg")
println("Moment necessary to impose the constraint at hinge node = $hingeBalanceM")

println("Finished sweptBackFFWTWingHingeAxisAndFoldConstraint.jl")