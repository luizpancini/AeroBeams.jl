using AeroBeams, LinearAlgebra

# Beam
L = 10
EA,GA,GJ,EI = 1e6,1e6,1e6,1e5
ρA,ρI = 1,1
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

# BCs 
τ = 0.5
θ₀ = -1.5
θ = t -> ifelse.(t.<=τ, θ₀*t/τ, θ₀)
p = t -> 4*tan.(θ(t)/4)
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
driver = create_BC(name="driver",beam=beam,node=1,types=["p2A"],values=[t->p(t)])

# Model
rotationDrivenArticulatedRobotArm = create_Model(name="rotationDrivenArticulatedRobotArm",beams=[beam],BCs=[support,driver])

# Time variables
tf = 2.2
Δt = 1e-2

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2, Δt=1e-6)

# Create and solve the problem
problem = create_DynamicProblem(model=rotationDrivenArticulatedRobotArm,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
u1_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[1] for i in 1:length(t)]
u3_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[3] for i in 1:length(t)]

println("Finished rotationDrivenArticulatedRobotArm.jl")