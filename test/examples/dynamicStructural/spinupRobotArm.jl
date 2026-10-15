using AeroBeams, ForwardDiff

# Select to promote rotation through basis A (does not work otherwise, perhaps because the rotation goes beyond the singularity of the rotation parameters)
rotateBasisA = true

# Rotation variables
τ = 15
θ = t -> ifelse(t<=τ, 6/τ*(t^2/2 + (τ/(2π))^2 * (cos(2π*t/τ)-1)), 6*t-45)
p = t -> 4 * tan(θ(t)/4) * [0,0,1]
ps = t -> scaled_rotation_parameters(p(t))
θdot = t -> ForwardDiff.derivative(θ,t)

# Beam
L = 10
EA,GA = 2.8e7, 1e7
GJ = EI = 1.4e4
ρA,ρI = 1.2, 6e-4
nElements = 5
beam = create_Beam(name="beam",length=L,nElements=nElements,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)])

# BCs
bc = rotateBasisA ? 
     create_BC(name="pin",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0)) : 
     create_BC(name="driver",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>t->ps(t)[1],:p2A=>t->ps(t)[2],:p3A=>t->ps(t)[3]))

# Angular velocity of basis A
ω_A = t -> rotateBasisA ? [0,0,-θdot(t)] : [0,0,0]

# Model
spinupRobotArm = create_Model(name="spinupRobotArm",beams=[beam],BCs=[bc],ω_A=ω_A)

# Time variables
tf = 20
Δt = 4e-2
timeVector = unique(vcat(collect(0:Δt/2:τ),collect(τ:Δt:tf)))

# Create and solve the problem
problem = create_DynamicProblem(model=spinupRobotArm,timeVector=timeVector)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][nElements].nodal.u2[1] for i in 1:length(t)]
u2_tip = [states[i][nElements].nodal.u2[2] for i in 1:length(t)]
θ3_root = [states[i][1].nodal.θ[1] for i in 1:length(t)]

# Rigid beam's tip displacements
u1_tip_rigid = L*(cos.(θ.(t)).-1)
u2_tip_rigid = L*sin.(θ.(t))

## Plots
# Set paths
relPath = "/test/outputs/figures/spinupRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,follow=:basisA,plotFrequency=5,showScale=false,plotLimits=([-L,L],[-L,L],[0,L]),save=true,savePath=string(relPath,"/spinupRobotArm_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
ms = 5
msw = 0
gr()

# Normalized tip u₁
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_1/L\$")
plot!(t,u1_tip/L, c=:black, lw=lw, label="Flexible")
scatter!(t[1:5:end],u1_tip_rigid[1:5:end]/L, c=:blue, ms=ms, msw=msw, label="Rigid")
display(plt1)
savefig(string(absPath,"/spinupRobotArm_u1.pdf"))

# Normalized tip u₂
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_2/L\$")
plot!(t,u2_tip/L, c=:black, lw=lw, label="Flexible")
scatter!(t[1:5:end],u2_tip_rigid[1:5:end]/L, c=:blue, ms=ms, msw=msw, label="Rigid")
display(plt2)
savefig(string(absPath,"/spinupRobotArm_u2.pdf"))

# Root rotation
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root \$\\theta/\\pi\$")
plot!(t,θ3_root/π, c=:black, lw=lw, label="Numerical")
scatter!(t[1:20:end],θ.(t[1:20:end])/π, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt3)
savefig(string(absPath,"/spinupRobotArm_theta.pdf"))

# Axial force
plot_time_outputs(problem,nodes=[(1,1)],elements=[1,nElements],nodalOutputs=["F1"],elementalOutputs=["F1"],save=true,saveFolder=string(relPath,"/"))

println("Finished spinupRobotArm.jl")