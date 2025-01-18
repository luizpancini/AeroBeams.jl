using AeroBeams, LinearAlgebra

# Spring stiffness range
kRange = vcat(5e-3,1e-2,1e-1,5e-1,1,10)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = -1
q = (x1,t) -> q₀
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# System solver
σ0 = 1
maxIter = 20
relTol = 1e-9
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Model
clampedHingedBeamSpringRange = create_Model(name="clampedHingedBeamSpringRange",beams=[beam],BCs=[clamp])

# Initialize outputs
u1 = Vector{Vector{Float64}}(undef,length(kRange))
u3 = Vector{Vector{Float64}}(undef,length(kRange))
p2 = Vector{Vector{Float64}}(undef,length(kRange))
F3 = Vector{Vector{Float64}}(undef,length(kRange))
M2 = Vector{Vector{Float64}}(undef,length(kRange))
springs = Vector{Spring}(undef,length(kRange))
hingeAngle = Vector{Float64}(undef,length(kRange))
springMoment = Vector{Float64}(undef,length(kRange))
problem = Vector{SteadyProblem}(undef,length(kRange))

# Loop spring stiffness
for (i,k) in enumerate(kRange)
    # Show progress
    println("Solving for k = $k")
    # Remove springs from beam
    remove_all_springs_from_beams!(beams=[beam])
    # Update spring stiffness on beam
    springs[i] = create_Spring(elementsIDs=[hingeNode-1,hingeNode],nodesSides=[1,2],kp=[0,k,0])
    add_spring_to_beams!(beams=[beam,beam],spring=springs[i])
    # Update model
    update_model!(clampedHingedBeamSpringRange)
    # Set initial guess solution as previous known solution
    x0 = i==1 ? zeros(0) : problem[i-1].x
    # Create and solve the problem
    problem[i] = create_SteadyProblem(model=clampedHingedBeamSpringRange,systemSolver=NR,x0=x0)
    solve!(problem[i])
    # Get outputs
    u1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
    u3[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    p2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].p_n1[2],problem[i].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
    F3[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].F_n1[3],problem[i].nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
    M2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].M_n1[2],problem[i].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
    hingeAngle[i] = rotation_angle_limited([0; p2[i][2*hingeNode-1]-p2[i][2*hingeNode-2]; 0])*180/pi
    springMoment[i] = norm(springs[i].Ms)
end

# Get nodal arclength positions
elemNodes = vcat([vcat(problem[1].model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem[1].model.r_n[n][1] for n in elemNodes]
r_n2 = [problem[1].model.r_n[n][2] for n in elemNodes]
r_n3 = [problem[1].model.r_n[n][3] for n in elemNodes]
x1 = vcat([vcat(problem[1].model.beams[1].elements[e].x1_n1,problem[1].model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

println("Finished clampedHingedBeamSpringRange.jl")