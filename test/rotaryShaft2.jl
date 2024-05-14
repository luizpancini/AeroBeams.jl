using AeroBeams, LinearAlgebra, ForwardDiff, Plots

# Rotation variables
Δθ = 1*π/180
ω = 1*π/180
T = 2*π/ω
θ = t -> Δθ*sin(ω*t)
θdot = t -> ForwardDiff.derivative(θ,t)
θddot = t -> ForwardDiff.derivative(θdot,t)
p = t -> 4*tan(θ(t)/4)

# Beam
L = 1
∞ = 1e12
nElem = 1
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1,ρIs=1e-4)])

# BCs
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])
journal = create_BC(name="journal",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
rotaryShaft = create_Model(name="rotaryShaft",beams=[beam],BCs=[support,journal],ω_A=t->[θdot(t);0;0])

# Time variables
cycles = 1
tf = cycles*T
Δt = T/1e2

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

# Analytical solution
out = t -> rotation_tensor_WM([p(t);0;0]) 
R = t -> out(t)[1]
HT = t -> tangent_operator_transpose_WM([p(t);0;0])
ΩAnalytical = t -> θdot(0)
ΩdotAnalytical = t -> θddot(0)
pdotVecAnalytical = t -> HT(t)'*([ΩAnalytical(t);0;0]-R(t)'*[θdot(t);0;0])
pdotAnalytical = t -> pdotVecAnalytical(t)[1]

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
plotRange = 1:5:length(t)
# Normalized rotation parameter
plt1 = plot(xlabel="\$t/T\$", ylabel="\$p_1/\\Delta\\theta\$ ")
plot!(tNorm,pNum/Δθ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[plotRange],p.(t[plotRange])/Δθ, c=:blue, ms=ms, label="Analytical")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft2_1.pdf"))
# Normalized rotation parameter rate 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{p}_1/\\Delta\\theta\\omega\$ [\$1\$/s]")
plot!(tNorm,pdotNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[plotRange],pdotAnalytical.(t[plotRange])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft2_2.pdf"))
# Normalized sectional angular velocity 
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1/\\Delta\\theta\\omega\$ [1/s]", ylims=[0,2])
plot!(tNorm,ΩNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[plotRange],ΩAnalytical.(t[plotRange])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft2_3.pdf"))
# Normalized sectional angular acceleration 
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1/\\Delta\\theta\\omega^2\$ [1/\$s^2\$]", ylims=[-1,1])
plot!(tNorm,ΩdotNum/(Δθ*ω^2), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[plotRange],ΩdotAnalytical.(t[plotRange])/(Δθ*ω^2), c=:blue, ms=ms, label="Analytical")
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/rotaryShaft2_4.pdf"))

println("Finished rotaryShaft2.jl")