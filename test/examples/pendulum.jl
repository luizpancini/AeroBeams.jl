using AeroBeams, LinearAlgebra

# Gravity
g = 9.80665

# Initial angle of release
θ₀ = 1/8 * π/2

# Beam
L,r = 1.0,0.01
A,J = π*r^2,π/2*r^4
I = J/2
E,G,ρ = 200e9,80e9,7.9e3
nElem = 5
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*I,E*I])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,2ρ*I,ρ*I,ρ*I])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,(π/2-θ₀),0])

# BCs
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])

# Model
pendulum = create_Model(name="pendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g])

# Time variables
ω = sqrt(3/2*g/L)
T = 2*π/ω
cycles = 1
tf = cycles*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=pendulum,finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][end].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

# Analytical solution (valid for small θ₀)
θ = θ₀*cos.(ω*t)
u1_tip_analytical = -L*(sin(θ₀) .- sin.(θ))
u3_tip_analytical = L*(cos(θ₀) .- cos.(θ))

println("Finished pendulum.jl")