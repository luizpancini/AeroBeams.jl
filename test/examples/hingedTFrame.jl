using AeroBeams, LinearAlgebra

# Choose to calculate linear or nonlinear solution (they differ more and more as the distributed load q increases)
linearSolution = true

# Beam frame
L = 1
EI = 200e9*5e-6
∞ = 1e12
stiffnessMatrix = diagm([∞,∞,∞,∞,EI,EI])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[stiffnessMatrix],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[true,false,true]])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[stiffnessMatrix],rotationParametrization="E321",p0=[π/2;0;0],hingedNodes=[1],hingedNodesDoF=[[true,false,true]],connectedBeams=[beam1],connectedNodesThis=[1],connectedNodesOther=[div(nElem,2)+1])

# BCs
q = 1e5
add_loads_to_beam!(beam2,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[q; 0; q]])
pin1 = create_BC(name="pin1",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
pin2 = create_BC(name="pin2",beam=beam1,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
pin3 = create_BC(name="pin3",beam=beam2,node=nElem+1,types=["u1A","u3A","p2A"],values=[0,0,0])

# Model
hingedTFrame = create_Model(name="hingedTFrame",beams=[beam1,beam2],BCs=[pin1,pin2,pin3])

# Setup nonlinear system solver
σ0 = 0.01
σstep = 0.01
atol = 1e-4
maxit = 50
NR = create_NewtonRaphson(absoluteTolerance=atol,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxit,displayStatus=true)

# Create and solve the problem
problem = create_SteadyProblem(model=hingedTFrame,systemSolver=NR,getLinearSolution=linearSolution)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1_beam1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
x1_beam2 = vcat([vcat(problem.model.beams[2].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u1_beam1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u3_beam1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3_beam1 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2_beam1 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
u1_beam2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in nElem+1:2*nElem]...)
u3_beam2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in nElem+1:2*nElem]...)

println("Finished hingedTFrame.jl")