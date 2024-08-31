using AeroBeams, LinearAlgebra, DelimitedFiles

# Beam frame
L = 120
A,Iy = 6,2
E = 7.2e6
ν = 0.0
G = E/(2*(1+ν))
∞ = 1e10
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-π/2;0])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
F = 18e3
elemForce = div(nElem,5)
support1 = create_BC(name="support1",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
support2 = create_BC(name="support2",beam=beam2,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
force = create_BC(name="force",beam=beam2,node=elemForce+1,types=["F3A"],values=[-F])

# Model
LeeFrameDeadLoad = create_Model(name="LeeFrameDeadLoad",beams=[beam1,beam2],BCs=[support1,support2,force],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=LeeFrameDeadLoad,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
u1_atForce = [problem.nodalStatesOverσ[i][nElem+elemForce].u_n2[1] for i in 1:length(σVector)]
u3_atForce = [problem.nodalStatesOverσ[i][nElem+elemForce].u_n2[3] for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm("test/referenceData/LeeFrameDeadLoad/u1.txt")
u3Ref = readdlm("test/referenceData/LeeFrameDeadLoad/u3.txt")

println("Finished LeeFrameDeadLoad.jl")