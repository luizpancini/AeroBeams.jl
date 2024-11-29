using AeroBeams, LinearAlgebra

# Beam
L = 1
∞ = 1e6
stiffnessMatrix = diagm([∞,∞,∞,∞,∞,∞])
nElem1 = 2
nElem2 = 1
beam1 = create_Beam(name="beam1",length=L,nElements=nElem1,S=[stiffnessMatrix])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem2,S=[stiffnessMatrix],rotationParametrization="E321",p0=[0;π/2;0])

# BCs - balanceLoads are the equivalent reactions at the pin
F = 1
roller = create_BC(name="roller",beam=beam1,node=1,types=["u3A"],values=zeros(1))
pin = create_BC(name="pin",beam=beam1,node=nElem1+1,types=["u1A","u2A","u3A","p1A","p3A"],values=zeros(5))
verticalForce = create_BC(name="verticalForce",beam=beam1,node=div(nElem1,2)+1,types=["F3A"],values=[F])
horizontalForce = create_BC(name="horizontalForce",beam=beam2,node=nElem2+1,types=["F1A"],values=[F])
balanceLoads = create_BC(name="balanceLoads",beam=beam2,node=1,types=["F1A","F3A"],values=zeros(2),toBeTrimmed=trues(2))

# Model
rightAngledFrameTrim = create_Model(name="rightAngledFrameTrim",beams=[beam1,beam2],BCs=[roller,pin,verticalForce,horizontalForce,balanceLoads])

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=100,displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=rightAngledFrameTrim,systemSolver=NR)
solve!(problem)

# Get solution 
balanceHorizontalForce = problem.x[end-1]*problem.model.forceScaling
balanceVerticalForce = problem.x[end]*problem.model.forceScaling 

# Compare to analytical solution 
balanceHorizontalForceAnalytical = -F
balanceVerticalForceAnalytical = -(F*L/2+F*L)/L

ϵ_rel_F1 = balanceHorizontalForce/balanceHorizontalForceAnalytical - 1
ϵ_rel_F3 = balanceVerticalForce/balanceVerticalForceAnalytical - 1

println("Relative errors:\nF1: $ϵ_rel_F1 \nF3: $ϵ_rel_F3")

println("Finished rightAngledFrameTrim.jl")