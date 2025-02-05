using AeroBeams, LinearAlgebra

# Beam rotation angle about a3-axis
α = π*30/180

# Beam rotation angle about a2'-axis
β = -π*0/180

# Hinge angle [rad]
θ = -π*45/180

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.25
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],rotationParametrization="E321",p0=[α;β;0],aeroSurface=surf)

# Hinge axis constraint
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=[0;1;0],pHValue=4*tan(θ/4))

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
clampedHingedBeamRotated = create_Model(name="clampedHingedBeamRotated",beams=[beam],BCs=[clamp,hingeForce,tipForce],hingeAxisConstraints=[hingeAxisConstraint])

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
p1_b = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1_b[1],problem.nodalStatesOverσ[end][e].p_n2_b[1]) for e in 1:nElem]...)
p2_b = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1_b[2],problem.nodalStatesOverσ[end][e].p_n2_b[2]) for e in 1:nElem]...)
p3_b = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1_b[3],problem.nodalStatesOverσ[end][e].p_n2_b[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Compute rotation, displacement and balance moment at the hinge node
p_bOutboard = [p1_b[2*hingeNode-1]; p2_b[2*hingeNode-1]; p3_b[2*hingeNode-1]]
p_bInboard = [p1_b[2*hingeNode-2]; p2_b[2*hingeNode-2]; p3_b[2*hingeNode-2]]
Δp_bHingeNode = rotation_between_WM(p_bInboard,p_bOutboard)
ypr = WM_to_ypr(Δp_bHingeNode)*180/π
u3HingeNode = u3[2*hingeNode-1]
u3Tip = u3[end]
hingeBalanceM = -problem.model.hingeAxisConstraints[1].balanceMoment

# Display results
println("Rotation at hinge node: yaw = $(ypr[1]) deg, pitch = $(ypr[2]) deg, roll = $(ypr[3])")
println("Displacement at hinge node = $u3HingeNode")
println("Displacement at tip = $u3Tip")
println("Moment necessary to impose the constraint at hinge node = $hingeBalanceM")

println("Finished clampedHingedBeamRotated.jl")