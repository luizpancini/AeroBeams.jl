using AeroBeams, LinearAlgebra

# Gravity
g = 9.80665

# Initial angle of release
θ₀ = π/2

# Beam 
L,b,H = 2,0.01,0.01
A,I = b*H,b*H^3/12
E,ρ = 200e9,7.8e3
∞ = 1e16
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*I,∞])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,2*ρ*I,ρ*I,ρ*I])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,(π/2-θ₀),0],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

# BCs
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])

# Model
doublePendulum = create_Model(name="doublePendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g])

# Time variables
tf = 1
Δt = tf/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=doublePendulum,finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[1] for i in 1:length(t)]
u3_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[3] for i in 1:length(t)]
u1_tip = [problem.nodalStatesOverTime[i][end].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

println("Finished doublePendulum.jl")