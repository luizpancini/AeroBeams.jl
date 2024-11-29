using AeroBeams, LinearAlgebra

# Beam
L,b,H = 1.0,0.1,0.1
E,ρ = 200e9,7.8e3
A,Iy,Iz, = b*H,b*H^3/12,H*b^3/12
Is = Iy+Iz
∞ = 1e12
nElements = 20
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*Iy,∞])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElements,S=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
F = 1e3
ω = 200
force = create_BC(name="force",beam=beam,node=nElements+1,types=["F3A"],values=[t->F*sin(ω*t)])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
tipSineLoadedCantilever = create_Model(name="tipSineLoadedCantilever",beams=[beam],BCs=[clamp,force])

# Time variables
T = 2*π/ω
tf = 20*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=tipSineLoadedCantilever,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]
F3_root = [problem.nodalStatesOverTime[i][1].F_n1[3] for i in 1:length(t)]
M2_root = [problem.nodalStatesOverTime[i][1].M_n1[2] for i in 1:length(t)]

println("Finished tipSineLoadedCantilever.jl")