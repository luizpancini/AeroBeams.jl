using AeroBeams, LinearAlgebra

# Option for linear solution
linear = false

# Spring stiffness range
kRange = vcat(1e-4,collect(1e-3:1e-3:1e-1),collect(1:1:10))

# TF for using previous solution as initial guess
usePreviousSol = trues(length(kRange))
usePreviousSol[end-9:end] .= false

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[false,true,false]])

# Initialize spring around hinge
spring = create_Spring(basis="A",elementsIDs=[hingeNode-1,hingeNode],nodesSides=[1,2],kOOPBending=kRange[1])
add_spring_to_beams!(beams=[beam,beam],spring=spring)

# BCs
q₀ = -1
q = (x1,t) -> q₀
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# System solver
σ0 = 1e-0
maxIter = 500
relTol = 1e-5
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Model
clampedHingedBeamSpringRange = create_Model(name="clampedHingedBeamSpringRange",beams=[beam],BCs=[clamp])

# Initialize outputs
u1 = Vector{Vector{Float64}}(undef,length(kRange))
u3 = Vector{Vector{Float64}}(undef,length(kRange))
p2 = Vector{Vector{Float64}}(undef,length(kRange))
F3 = Vector{Vector{Float64}}(undef,length(kRange))
M2 = Vector{Vector{Float64}}(undef,length(kRange))
problem = Vector{SteadyProblem}(undef,length(kRange))

# Loop spring stiffness
for (i,k) in enumerate(kRange)
    # Show progress
    println("Solving for k = $k")
    # Remove springs from beam
    clampedHingedBeamSpringRange.beams[1].springs = Vector{Spring}()
    # Update spring stiffness on beam
    spring = create_Spring(basis="A",elementsIDs=[hingeNode-1,hingeNode],nodesSides=[1,2],kOOPBending=k)
    add_spring_to_beams!(beams=[beam,beam],spring=spring)
    # Update model
    update_model!(clampedHingedBeamSpringRange)
    # Set initial guess solution as previous known solution
    if i > 1 && usePreviousSol[i] && problem[i-1].systemSolver.convergedFinalSolution
        x0 = problem[i-1].x
    else
        x0 = zeros(0)
    end
    # Create and solve the problem
    problem[i] = create_SteadyProblem(model=clampedHingedBeamSpringRange,getLinearSolution=linear,systemSolver=NR,x0=x0)
    solve!(problem[i])
    # TF for converged solution
    converged = problem[i].systemSolver.convergedFinalSolution
    # Get outputs
    u1[i] = converged ? vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...) : fill(NaN,2*nElem)
    u3[i] = converged ? vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...) : fill(NaN,2*nElem)
    p2[i] = converged ? vcat([vcat(scaled_rotation_parameters(problem[i].nodalStatesOverσ[end][e].p_n1)[2],scaled_rotation_parameters(problem[i].nodalStatesOverσ[end][e].p_n2)[2]) for e in 1:nElem]...) : fill(NaN,2*nElem)
    F3[i] = converged ? vcat([vcat(problem[i].nodalStatesOverσ[end][e].F_n1[3],problem[i].nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...) : fill(NaN,2*nElem)
    M2[i] = converged ? vcat([vcat(problem[i].nodalStatesOverσ[end][e].M_n1[2],problem[i].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...) : fill(NaN,2*nElem)
end

# Get nodal arclength positions
elemNodes = vcat([vcat(problem[1].model.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [problem[1].model.r_n[n][1] for n in elemNodes]
r_n2 = [problem[1].model.r_n[n][2] for n in elemNodes]
r_n3 = [problem[1].model.r_n[n][3] for n in elemNodes]
x1 = vcat([vcat(problem[1].model.beams[1].elements[e].x1_n1,problem[1].model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

println("Finished clampedHingedBeamSpringRange.jl")