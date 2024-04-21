using AeroBeams, LinearAlgebra, Plots

# Rotation variables
θ = π/3
ω = 30
T = 2*π/ω
θ₁ = t -> θ*sin.(ω*t)
θdot₁ = t -> θ*ω*cos.(ω*t)
θddot₁ = t -> -θ*ω^2*sin.(ω*t)
p₁ = t -> 4*tan.(θ₁(t)/4)
pdot₁ = t -> sec.(θ₁(t)/4).^2 .* θdot₁(t)
pddot₁ = t -> sec.(θ₁(t)/4).^2 .* (1/8*p₁(t).*θdot₁(t).^2 .+ θddot₁(t))

# Beam
L,r = 1.0,0.05
A,J = π*r^2,π/2*r^4
Iy,Iz,Is = J/2,J/2,J
E = 200e9
G = E/(2*(1+0.3))
ρ = 7.9e3
nElem = 5
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*Iy,E*Iz])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],pdot0_of_x1=x1->[pdot₁(0); 0.0; 0.0])

# BCs
driver = create_BC(name="driver",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t -> p₁(t),0,0])
journal = create_BC(name="journal",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
rotaryShaft = create_Model(name="rotaryShaft",beams=[beam],BCs=[driver,journal])

# Time variables
cycles = 1
tf = cycles*T
Δt = 5e-4

# Create and solve the problem
problem = create_DynamicProblem(model=rotaryShaft,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
e = 1
p₁_num = [problem.elementalStatesOverTime[i][e].p[1] for i in 1:length(tNorm)]
pdot₁_num = [problem.elementalStatesRatesOverTime[i][e].pdot[1] for i in 1:length(tNorm)]
pddot₁_num = [problem.elementalStatesRatesOverTime[i][e].pddot[1] for i in 1:length(tNorm)]
Ω₁_num = [problem.elementalStatesOverTime[i][e].Ω[1] for i in 1:length(tNorm)]
Ωdot₁_num = [problem.elementalStatesRatesOverTime[i][e].Ωdot[1] for i in 1:length(tNorm)]
M₁_num = [problem.elementalStatesOverTime[i][1].M[1] for i in 1:length(tNorm)]

# Plots
# ------------------------------------------------------------------------------
# Normalized rotation parameter
plt1 = Plots.plot()
Plots.plot!(tNorm,p₁_num/θ, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$p_1/\\theta\$ ", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],p₁(t[1:20:end])/θ, c=:blue, markersize=3, label="Analytical", show=true)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_1.pdf"))
# Normalized rotation parameter rate 
plt2 = Plots.plot()
Plots.plot!(tNorm,pdot₁_num/(θ*ω), c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\dot{p}_1/\\theta\\omega\$ [\$1\$/s]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],pdot₁(t[1:20:end])/(θ*ω), c=:blue, markersize=3, label="Analytical", show=true)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_2.pdf"))
# Normalized rotation parameter second rate
plt3 = Plots.plot()
Plots.plot!(tNorm,pddot₁_num/(θ*ω^2), c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\ddot{p}_1/\\theta\\omega^2\$ [\$1\$/\$s^2\$]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],pddot₁(t[1:20:end])/(θ*ω^2), c=:blue, markersize=3, label="Analytical", show=true)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_3.pdf"))
# Normalized angular velocity 
plt4 = Plots.plot()
Plots.plot!(tNorm,Ω₁_num/(θ*ω), c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\Omega_1/\\theta\\omega\$ [1/s]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],θdot₁(t[1:20:end])/(θ*ω), c=:blue, markersize=3, label="Analytical", show=true)
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_4.pdf"))
# Normalized angular acceleration 
plt5 = Plots.plot()
Plots.plot!(tNorm,Ωdot₁_num/(θ*ω^2), c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1/\\theta\\omega^2\$ [1/\$s^2\$]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],θddot₁(t[1:20:end])/(θ*ω^2), c=:blue, markersize=3, label="Analytical", show=true)
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_5.pdf"))
# Normalized driving torque
plt6 = Plots.plot()
Plots.plot!(tNorm,M₁_num, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$M_1^*\$ [N.m]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],-θddot₁(t[1:20:end])*(ρ*Is), c=:blue, markersize=3, label="Analytical", show=true)
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft_6.pdf"))

println("Finished rotaryShaft.jl")