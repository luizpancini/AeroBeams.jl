using AeroBeams

# Gravity
g = standard_atmosphere(0).g

# Initial angle of release [rad]
θ₀ = 1/8 * π/2

# Beam
L,r = 1.0, 0.01
A,J = π*r^2, π/2*r^4
I = J/2
E,G,ρ = 200e9, 80e9, 7.9e3
nElem = 5
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,GAy=G*A,GAz=G*A,GJ=G*J,EIy=E*I,EIz=E*I)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*I,ρIz=ρ*I)],p0=transform_to_WM([0,(π/2-θ₀),0],from=:E321))

# BCs
support = create_BC(name="support",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))

# Model
pendulum = create_Model(name="pendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g])

# Time variables
ω = sqrt(3/2*g/L)
T = 2π/ω
cycles = 1
tf = cycles * T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=pendulum,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][end].nodal.u2[1] for i in 1:length(t)]
u3_tip = [states[i][end].nodal.u2[3] for i in 1:length(t)]

# Analytical solution (valid for small θ₀)
θ = θ₀*cos.(ω*t)
u1_tip_analytical = -L*(sin(θ₀) .- sin.(θ))
u3_tip_analytical = L*(cos(θ₀) .- cos.(θ))

## Plots
# Print warning
if abs(θ₀) > π/8
    @warn "Initial angle of release is large, analytical comparison is not fair"
end

# Set paths
relPath = "/test/outputs/figures/pendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=([-L,L],[-L,L],[-L,0]),save=true,savePath=string(relPath,"/pendulum_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
ms = 4
msw = 0
gr()

# Normalized tip u1 displacement
plt1 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_1/L\$")
plot!(t/T,u1_tip/L, c=:black, lw=lw, label="Numerical")
scatter!(t[1:2:end]/T,u1_tip_analytical[1:2:end]/L, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/pendulum_u1.pdf"))

# Normalized tip u3 displacement
plt2 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_3/L\$")
plot!(t/T,u3_tip/L, c=:black, lw=lw, label="Numerical")
scatter!(t[1:2:end]/T,u3_tip_analytical[1:2:end]/L, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/pendulum_u3.pdf"))

println("Finished pendulum.jl")