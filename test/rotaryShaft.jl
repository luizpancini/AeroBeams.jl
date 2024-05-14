using AeroBeams, LinearAlgebra, ForwardDiff, Plots

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
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
e = 1
pNum = [problem.elementalStatesOverTime[i][e].p[1] for i in 1:length(tNorm)]
pdotNum = [problem.elementalStatesRatesOverTime[i][e].pdot[1] for i in 1:length(tNorm)]
ΩNum = [problem.elementalStatesOverTime[i][e].Ω[1] for i in 1:length(tNorm)]
ΩdotNum = [problem.elementalStatesRatesOverTime[i][e].Ωdot[1] for i in 1:length(tNorm)]
MNum = [problem.elementalStatesOverTime[i][1].M[1] for i in 1:length(tNorm)]

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Normalized rotation parameter
plt1 = plot(xlabel="\$t/T\$", ylabel="\$p_1/\\Delta\\theta\$ ")
plot!(tNorm,pNum/Δθ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],p.(t[1:20:end])/Δθ, c=:blue, ms=ms, label="Analytical")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_1.pdf"))
# Normalized rotation parameter rate 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{p}_1/\\Delta\\theta\\omega\$ [\$1\$/s]")
plot!(tNorm,pdotNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],pdot.(t[1:20:end])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_2.pdf"))
# Normalized sectional angular velocity 
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1/\\Delta\\theta\\omega\$ [1/s]")
plot!(tNorm,ΩNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],θdot.(t[1:20:end])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_4.pdf"))
# Normalized sectional angular acceleration 
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1/\\Delta\\theta\\omega^2\$ [1/\$s^2\$]")
plot!(tNorm,ΩdotNum/(Δθ*ω^2), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],θddot.(t[1:20:end])/(Δθ*ω^2), c=:blue, ms=ms, label="Analytical")
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_5.pdf"))
# Normalized driving torque
plt6 = plot(xlabel="\$t/T\$", ylabel="\$M_1^*\$ [N.m]")
plot!(tNorm,MNum, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],-θddot.(t[1:20:end])*(ρ*Is), c=:blue, ms=ms, label="Analytical")
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_6.pdf"))

println("Finished rotaryShaft.jl")