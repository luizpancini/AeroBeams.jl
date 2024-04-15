using AeroBeams, LinearAlgebra, Plots

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,1e3,1e3
ρA,ρI = 1,10
nElem = 40
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
M₀ = -80
τ = 2.5-0.1
M = t -> ifelse.(t.<=τ, M₀, 0)
hinge = create_BC(name="hinge",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
moment = create_BC(name="moment",beam=beam,node=1,types=["M2A"],values=[t->M(t)])

# Model
momentDrivenRobotArm = Model(name="momentDrivenRobotArm",beams=[beam],BCs=[hinge,moment])

# Time variables
tf = 15
Δt = 1e-1

# Create and solve the problem
problem = DynamicProblem(model=momentDrivenRobotArm,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u₁_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u₃_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
# Nomalized tip u₁
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,u₁_tip/L, linewidth=2, label="\$u_1/L\$")
plot!(t,u₃_tip/L, linewidth=2, label="\$u_3/L\$")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/momentDrivenRobotArm.pdf"))

println("Finished momentDrivenRobotArm.jl")