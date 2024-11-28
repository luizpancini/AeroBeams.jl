using AeroBeams, LinearAlgebra

# Beam
R,θ = 100,π/4
L = R*θ
EA,GA,GJ,EI = 1e7,5e6,833_333,833_333
ρA,ρI = 1,10
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],k=[0;0;1/R])

# BCs
F₀ = 100
F = t -> F₀*t
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[(t)->F(t)])

# Model
curvedCantileverDynamicFollower = create_Model(name="curvedCantileverDynamicFollower",beams=[beam],BCs=[clamp,tipForce])

# Time variables
tf = 7.5
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=curvedCantileverDynamicFollower,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u2_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[2] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
F1_root = [problem.nodalStatesOverTime[i][1].F_n1[1] for i in 1:length(t)]
F2_root = [problem.nodalStatesOverTime[i][1].F_n1[2] for i in 1:length(t)]
F3_root = [problem.nodalStatesOverTime[i][1].F_n1[3] for i in 1:length(t)]
M1_root = [problem.nodalStatesOverTime[i][1].M_n1[1] for i in 1:length(t)]
M2_root = [problem.nodalStatesOverTime[i][1].M_n1[2] for i in 1:length(t)]
M3_root = [problem.nodalStatesOverTime[i][1].M_n1[3] for i in 1:length(t)]

println("Finished curvedCantileverDynamicFollower.jl")
