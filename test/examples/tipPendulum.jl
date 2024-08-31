using AeroBeams, LinearAlgebra

# Gravity
g = 9.81

# Initial angle of release
θ₀ = 1/8 * π/2

# Beam
L,r = 1,0.01
A,J = π*r^2,π/2*r^4
I,Is = J/2,J
E,G,ρ = 200e9,80e9,7.9e3
nElem = 10
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*I,E*I])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*I,ρ*I])
beamMass = ρ*A*L
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,(π/2-θ₀),0])

# Pendulum's tip mass
massRatio = 10
tipMass = PointInertia(elementID=nElem,η=[L/nElem/2;0;0],mass=massRatio*beamMass)
add_point_inertias_to_beam!(beam,inertias=[tipMass])

# BCs
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])

# Model
tipPendulum = create_Model(name="tipPendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g])

# Time variables
ω = sqrt(g/L)
T = 2*π/ω
cycles = 2
tf = cycles*T
Δt = T/100

# Create and solve the problem (skip initial states update, otherwise it will set the pendulum to an equilibrium position)
problem = create_DynamicProblem(model=tipPendulum,finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][end].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

# Analytical solution (valid for small θ₀)
θ = θ₀*cos.(ω*t)
u1_tip_analytical = -L*(sin(θ₀) .- sin.(θ))
u3_tip_analytical = L*(cos(θ₀) .- cos.(θ))

println("Finished tipPendulum.jl")