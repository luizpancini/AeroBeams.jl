using AeroBeams, LinearAlgebra, ForwardDiff

# Rotation variables
Δθ = π/3
ω = 30
T = 2*π/ω
θ = t -> Δθ*sin(ω*t)
θdot = t -> ForwardDiff.derivative(θ,t)
θddot = t -> ForwardDiff.derivative(θdot,t)
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)
pddot = t -> ForwardDiff.derivative(pdot,t)

# Beam
L,r = 1.0,0.05
A,J = π*r^2,π/2*r^4
Iy,Iz,Is = J/2,J/2,J
E = 200e9
G = E/(2*(1+0.3))
ρ = 7.9e3
nElem = 10
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*Iy,E*Iz])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],pdot0_of_x1=x1->[pdot(0); 0.0; 0.0])

# BCs
driver = create_BC(name="driver",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t->p(t),0,0])
journal = create_BC(name="journal",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
rotaryShaft = create_Model(name="rotaryShaft",beams=[beam],BCs=[driver,journal])

# Time variables
cycles = 1
tf = cycles*T
Δt = T/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=rotaryShaft,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
e = 1
pNum = [problem.elementalStatesOverTime[i][e].p[1] for i in 1:length(tNorm)]
pdotNum = [problem.elementalStatesRatesOverTime[i][e].pdot[1] for i in 1:length(tNorm)]
ΩNum = [problem.elementalStatesOverTime[i][e].Ω[1] for i in 1:length(tNorm)]
ΩdotNum = [problem.elementalStatesRatesOverTime[i][e].Ωdot[1] for i in 1:length(tNorm)]
MNum = [problem.elementalStatesOverTime[i][1].M[1] for i in 1:length(tNorm)]

println("Finished rotaryShaft.jl")