using AeroBeams, LinearAlgebra

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,5e2,5e2
ρA,ρI = 1,10
θ₀ = atan(4/3)
nElem = 10
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0;θ₀;0])

# BCs - stable up to 18 s with M₀ = 100
M₀ = 200
τ = 2.5
M2 = t -> ifelse.(t.<=τ, M₀*t/τ, ifelse.(t.<=2*τ, 2*M₀*(1-t/(2*τ)), 0))
M3 = t -> M2(t)/2
F1 = t -> M2(t)/10
loads = create_BC(name="loads",beam=beam,node=nElem+1,types=["F1A","M2A","M3A"],values=[t->F1(t),t->M2(t),t->M3(t)])

# Model
flyingSpaghetti3D = create_Model(name="flyingSpaghetti3D",beams=[beam],BCs=[loads])

# Time variables
tf = 3.5
Δt = 1e-3

# Create and solve the problem
problem = create_DynamicProblem(model=flyingSpaghetti3D,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.savedTimeVector
u1_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in eachindex(t)]
u2_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[2] for i in eachindex(t)]
u3_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in eachindex(t)]
θ_tip = [problem.nodalStatesOverTime[i][nElem].θ_n2 for i in eachindex(t)]

println("Finished flyingSpaghetti3D.jl")