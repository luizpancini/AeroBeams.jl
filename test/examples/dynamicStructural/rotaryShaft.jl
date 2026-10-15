using AeroBeams, ForwardDiff

# Flag to apply a transition function to pitch profile
applyTransitionFunction = true

# Rotation variables
Δθ = π/3
ω = 30
τ = 2π/ω
θ = applyTransitionFunction ? 
    (t -> Δθ * sin(ω*t) * (t/τ)^10 / (1 + (t/τ)^10)) :
    (t -> Δθ * sin(ω*t))
θdot = t -> ForwardDiff.derivative(θ,t)
θddot = t -> ForwardDiff.derivative(θdot,t)
p = t -> 4 * tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)
pddot = t -> ForwardDiff.derivative(pdot,t)

# Beam
L,r = 1.0, 0.05
A,J = π*r^2, π/2*r^4
Iy = Iz = J/2
Is = Iy + Iz
E = 200e9
G = E/(2*(1+0.3))
ρ = 7.9e3
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,GAy=G*A,GAz=G*A,GJ=G*J,EIy=E*Iy,EIz=E*Iz)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)])

# Initial conditions
add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:pdot0_of_x1 => x1->[pdot(0), 0, 0]))

# BCs
driver = create_BC(name="driver",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>t->p(t),:p2A=>0,:p3A=>0))
journal = create_BC(name="journal",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))

# Model
rotaryShaft = create_Model(name="rotaryShaft",beams=[beam],BCs=[driver,journal])

# Time variables
T = 2π/ω
cycles = 3
tf = cycles * T
Δt = T/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=rotaryShaft,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
tNorm = t/T
pNum = [states[i][1].midpoint.p[1] for i in 1:length(t)]
ΩNum = [states[i][1].midpoint.Ω[1] for i in 1:length(t)]
MNum = [states[i][1].midpoint.M[1] for i in 1:length(t)]
pdotNum = [states[i][1].rates.pdot[1] for i in 1:length(t)]
ΩdotNum = [states[i][1].rates.Ωdot[1] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/rotaryShaft"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 3
gr()

# Normalized rotation parameter
plt1 = plot(xlabel="\$t/T\$", ylabel="\$p_1/\\Delta\\theta\$ ")
plot!(tNorm,pNum/Δθ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],p.(t[1:20:end])/Δθ, c=:blue, ms=ms, label="Analytical")
display(plt1)
savefig(string(absPath,"/rotaryShaft_p.pdf"))

# Normalized rotation parameter rate 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{p}_1/\\Delta\\theta\\omega\$ [\$1\$/s]")
plot!(tNorm,pdotNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],pdot.(t[1:20:end])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt2)
savefig(string(absPath,"/rotaryShaft_pdot.pdf"))

# Normalized sectional angular velocity 
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1/\\Delta\\theta\\omega\$ [1/s]")
plot!(tNorm,ΩNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],θdot.(t[1:20:end])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt4)
savefig(string(absPath,"/rotaryShaft_angVel.pdf"))

# Normalized sectional angular acceleration 
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1/\\Delta\\theta\\omega^2\$ [1/\$s^2\$]")
plot!(tNorm,ΩdotNum/(Δθ*ω^2), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],θddot.(t[1:20:end])/(Δθ*ω^2), c=:blue, ms=ms, label="Analytical")
display(plt5)
savefig(string(absPath,"/rotaryShaft_angAcc.pdf"))

# Normalized driving torque
plt6 = plot(xlabel="\$t/T\$", ylabel="\$M_1^*\$ [N.m]")
plot!(tNorm,MNum, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],-θddot.(t[1:20:end])*(ρ*Is), c=:blue, ms=ms, label="Analytical")
display(plt6)
savefig(string(absPath,"/rotaryShaft_torque.pdf"))

println("Finished rotaryShaft.jl")