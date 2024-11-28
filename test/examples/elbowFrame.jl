using AeroBeams, LinearAlgebra, DelimitedFiles

# Beam frame
L = 10
EA,GAy,GAz,GJ,EIy,EIz = 1e6,1e6,1e6,1e3,1e3,1e3
ρA,ρI = 1,10
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[π/2;0;0])

# BCs
F₀ = 50
F = t -> ifelse.(t.<=1, F₀*t/1, ifelse.(t.<=2, F₀*(2-t), 0.0))
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
elbowForce = create_BC(name="elbowForce",beam=beam2,node=1,types=["F3A"],values=[t->F(t)])

# Model
elbowFrame = create_Model(name="elbowFrame",beams=[beam1,beam2],BCs=[clamp,elbowForce])

# Time variables
tf = 30
Δt = 5e-2

# Create and solve the problem
problem = create_DynamicProblem(model=elbowFrame,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
t = problem.timeVector
u3_elbow = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

# Load reference solution
u3ElbowRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "elbowFrame", "u3_elbow.txt"))
u3TipRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "elbowFrame", "u3_tip.txt"))

println("Finished elbowFrame.jl")