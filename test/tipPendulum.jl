using AeroBeams, LinearAlgebra, Plots

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
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][end].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

# Analytical solution (valid for small θ₀)
θ = θ₀*cos.(ω*t)
u1_tip_analytical = -L*(sin(θ₀) .- sin.(θ))
u3_tip_analytical = L*(cos(θ₀) .- cos.(θ))

# Plots
# ------------------------------------------------------------------------------
# Normalized tip u1 displacement
plt1 = plot()
plot!(t/T,u1_tip/L, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="Tip \$u_1/L\$ ", label="Numerical")
scatter!(t[1:2:end]/T,u1_tip_analytical[1:2:end]/L, c=:blue, markersize=3, label="Analytical")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/tipPendulum_1.pdf"))
# Normalized tip u3 displacement
plt2 = plot()
plot!(t/T,u3_tip/L, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="Tip \$u_3/L\$ ", label="Numerical")
scatter!(t[1:2:end]/T,u3_tip_analytical[1:2:end]/L, c=:blue, markersize=3, label="Analytical")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/tipPendulum_2.pdf"))
# Gif of tip displacements on the plane
plt3 = plot()
r1_tip = problem.model.elements[end].r_n2[1]
r3_tip = problem.model.elements[end].r_n2[3]
anim = @animate for i=1:length(t)
    x = [(r1_tip + u1_tip[i])/L]
    y = [(r3_tip + u3_tip[i])/L]
    scatter(x, y, c=:black, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$", xlims=(-1, 1), ylims=(-1, 1), label=false)
end
gif(anim, string(pwd(),"/test/outputs/figures/tipPendulum.gif"), fps = 30)

if abs(θ₀) > π/8
    println("Initial angle of release is large, analytical comparison is not valid")
end

println("Finished tipPendulum.jl")