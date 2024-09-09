# # Dynamic analysis of an articulated beam
# This example illustrates how to set up a dynamic analysis, using the articulated beam in free flight proposed by [Simo and Vu-Quoc](https://doi.org/10.1115/1.3171871):

#md # ![](assets/articulatedbeam.png)
#md # *Articulated beam: definition and motion*

# ### Problem setup
# The set up of this problem involves the definition of the articulated links (beams).
using AeroBeams, LinearAlgebra

## Beam
L = 10
EA,GA,GJ,EI = 1e6,1e6,1e6,1e4
ρA,ρI1,ρI2 = 1,1,10
θ₀ = atan(4/3)
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix1 = diagm([ρA,ρA,ρA,2*ρI1,ρI1,ρI1])
inertiaMatrix2 = diagm([ρA,ρA,ρA,2*ρI2,ρI2,ρI2])
inertiaMatrices = vcat([inertiaMatrix2 for _ in 1:div(nElem,2)],[inertiaMatrix1 for _ in 1:div(nElem,2)])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=inertiaMatrices,rotationParametrization="E321",p0=[0;θ₀;0],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

## BCs 
M₀ = 160
τ = 0.5
M2 = t -> ifelse.(t.<=τ, M₀, 0)
F1 = t -> M2(t)/4
forces = create_BC(name="forces",beam=beam,node=nElem+1,types=["F1A","M2A"],values=[t->F1(t),t->M2(t)])

## Model
flyingScissors = create_Model(name="flyingScissors",beams=[beam],BCs=[forces])

## Time variables
tf = 5
Δt = 5e-2

## Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2, Δt=Δt/10)

## Create and solve the problem
problem = create_DynamicProblem(model=flyingScissors,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)

## Unpack numerical solution
t = problem.timeVector
u1_tipA = [problem.nodalStatesOverTime[i][1].u_n2[1] for i in 1:length(t)]
u3_tipA = [problem.nodalStatesOverTime[i][1].u_n2[3] for i in 1:length(t)]
u1_tipB = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u3_tipB = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
u1_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[1] for i in 1:length(t)]
u3_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[3] for i in 1:length(t)]

println("Finished flyingScissors.jl") #src