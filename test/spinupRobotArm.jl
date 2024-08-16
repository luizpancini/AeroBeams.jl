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

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 5
relPath = "/test/outputs/figures/spinupRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,refBasis="I",plotFrequency=5,showScale=false,timeStampPos=[0.1;-0.05;0],plotLimits=[(-L,L),(-L,L),(0,L)],save=true,savePath=string(relPath,"/spinupRobotArm_deformation.gif"),displayProgress=true)
# Normalized tip u₁
gr()
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_1/L\$")
plot!(t,u₁_tip/L, c=:black, lw=lw, label="Flexible")
scatter!(t[1:5:end],u₁_tip_rigid[1:5:end]/L, c=:blue, ms=ms, msw=0, label="Rigid")
display(plt1)
savefig(string(absPath,"/spinupRobotArm_u1.pdf"))
# Normalized tip u₂
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_2/L\$")
plot!(t,u₂_tip/L, c=:black, lw=lw, label="Flexible")
scatter!(t[1:5:end],u₂_tip_rigid[1:5:end]/L, c=:blue, ms=ms, msw=0, label="Rigid")
display(plt2)
savefig(string(absPath,"/spinupRobotArm_u2.pdf"))
# Root rotation
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root \$\\theta/\\pi\$")
plot!(t,θ₃_root/π, c=:black, lw=lw, label="Numerical")
scatter!(t[1:20:end],θ(t[1:20:end])/π, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt3)
savefig(string(absPath,"/spinupRobotArm_theta.pdf"))
# Axial force
plot_time_outputs(problem,nodes=[(1,1)],elements=[1,nElements],nodalOutputs=["F1"],elementalOutputs=["F1"],save=true,saveFolder=string(relPath,"/"))

println("Finished spinupRobotArm.jl")