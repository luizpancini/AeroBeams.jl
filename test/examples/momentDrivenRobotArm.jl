using AeroBeams, LinearAlgebra

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,1e3,1e3
ρA,ρI = 1,10
nElem = 40
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
M₀ = -80
τ = 2.5-0.1
M = t -> ifelse.(t.<=τ, M₀, 0)
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
moment = create_BC(name="moment",beam=beam,node=1,types=["M2A"],values=[t->M(t)])

# Model
momentDrivenRobotArm = create_Model(name="momentDrivenRobotArm",beams=[beam],BCs=[pin,moment])

# Time variables
tf = 15
Δt = 1e-1

# Create and solve the problem
problem = create_DynamicProblem(model=momentDrivenRobotArm,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]

println("Finished momentDrivenRobotArm.jl")