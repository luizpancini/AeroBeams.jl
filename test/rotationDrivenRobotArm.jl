using AeroBeams, LinearAlgebra, Plots

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,1e3,1e3
ρA,ρI = 1,10
nElem = 40
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
τ = 5
θ₀ = -1.5
θ = t -> ifelse.(t.<=τ, θ₀*t/τ, θ₀)
p = t -> 4*tan.(θ(t)/4)
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
driver = create_BC(name="driver",beam=beam,node=1,types=["p2A"],values=[t->p(t)])

# Model
rotationDrivenRobotArm = create_Model(name="rotationDrivenRobotArm",beams=[beam],BCs=[support,driver])

# Time variables
tf = 9
Δt = 1e-1

# Create and solve the problem
problem = create_DynamicProblem(model=rotationDrivenRobotArm,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u₁_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u₃_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
lw = 2
relPath = "/test/outputs/figures/rotationDrivenRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=[(-L,L),(0,L),(0,L)],save=true,savePath=string(relPath,"/rotationDrivenRobotArm_deformation.gif"),displayProgress=true)
# Nomalized tip displacements
gr()
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,u₁_tip/L, lw=lw, label="\$u_1/L\$")
plot!(t,u₃_tip/L, lw=lw, label="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/rotationDrivenRobotArm_disp.pdf"))

println("Finished rotationDrivenRobotArm.jl")