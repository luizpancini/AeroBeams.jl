using AeroBeams, LinearAlgebra

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,5e2,5e2
ρA,ρI = 1,10
θ₀ = atan(4/3)
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0;θ₀;0])

# BCs 
M₀ = 80
τ = 2.5
M2 = t -> ifelse.(t.<=τ, M₀, 0)
F1 = t -> M2(t)/10
loads = create_BC(name="loads",beam=beam,node=nElem+1,types=["F1A","M2A"],values=[t->F1(t),t->M2(t)])

# Model
flyingFlexibleBeam2D = create_Model(name="flyingFlexibleBeam2D",beams=[beam],BCs=[loads])

# Time variables
tf = 13
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=flyingFlexibleBeam2D,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]

println("Finished flyingFlexibleBeam2D.jl")