using AeroBeams, LinearAlgebra, Plots

# Select whether to promote rotation through basis A
rotateBasisA = true

# Rotation variables
τ = 15.0
θ = t -> ifelse.(t.<=τ, 6/15*(t.^2/2 .+ (15/(2*π))^2 * (cos.(2*π*t/τ).-1)), 6*t.-45)
θdot = t -> ifelse.(t.<=τ, 2*t/5 .- 3/π*sin.(2*pi*t/τ), 6)
p = t -> 4*tan.(θ(t)/4)
pdot = t -> sec.(θ(t)/4).^2 .* θdot(t)

# Beam
L = 10
EA,GA,GJ,EI = 2.8e7,1e7,1.4e4,1.4e4
ρA,ρI = 1.2,6e-4
nElements = 3
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElements,C=[stiffnessMatrix],I=[inertiaMatrix])
if !rotateBasisA
    add_initial_displacements_and_velocities_to_beam!(beam,conditionTypes=["pdot0_of_x1"],conditionFuns=[(x1)->[0; 0; pdot(0)]])
end

# BCs
driver = create_BC(name="driver",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,t->p(t)])
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A"],values=[0,0,0,0,0])

# Model
if rotateBasisA
    spinupRobotArm = create_Model(name="spinupRobotArm",beams=[beam],BCs=[pin],ω_A=t->[0;0;-θdot(t)])
else
    spinupRobotArm = create_Model(name="spinupRobotArm",beams=[beam],BCs=[driver])
end

# Time variables
tf = 20
Δt = 4e-2
time = unique(vcat(collect(0:Δt/2:tf/2),collect(tf/2:Δt:tf)))

# Create and solve the problem
problem = create_DynamicProblem(model=spinupRobotArm,timeVector=time)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u₁_tip = [problem.nodalStatesOverTime[i][nElements].u_n2[1] for i in 1:length(t)]
u₂_tip = [problem.nodalStatesOverTime[i][nElements].u_n2[2] for i in 1:length(t)]
θ₃_root = [problem.nodalStatesOverTime[i][1].θ_n1 for i in 1:length(t)]

# Rigid beam's tip displacements
u₁_tip_rigid = L*(cos.(θ(t)).-1)
u₂_tip_rigid = L*sin.(θ(t))

# Plots: flexible vs. rigid
# ------------------------------------------------------------------------------
# Normalized tip u₁
plt1 = plot()
plot!(t,u₁_tip/L, c=:black, linewidth=2, xlabel="\$t\$ [s]", ylabel="Tip \$u_1/L\$ ", label="Flexible", show=true)
scatter!(t[1:5:end],u₁_tip_rigid[1:5:end]/L, c=:blue, markersize=3, label="Rigid", show=true)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/spinupRobotArm_1.pdf"))
# Normalized tip u₂
plt2 = plot()
plot!(t,u₂_tip/L, c=:black, linewidth=2, xlabel="\$t\$ [s]", ylabel="Tip \$u_2/L\$", label="Flexible", show=true)
scatter!(t[1:5:end],u₂_tip_rigid[1:5:end]/L, c=:blue, markersize=3, label="Rigid", show=true)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/spinupRobotArm_2.pdf"))
# Root rotation
plt3 = plot()
plot!(t,θ₃_root/π, c=:black, linewidth=2, xlabel="\$t\$ [s]", ylabel="Root \$\\theta/\\pi\$", label="Numerical", show=true)
scatter!(t[1:20:end],θ(t[1:20:end])/π, c=:blue, markersize=3, label="Analytical", show=true)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/spinupRobotArm_3.pdf"))

println("Finished spinupRobotArm.jl")