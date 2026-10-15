using AeroBeams

# Beam
L = 10
EA = GA = 1e4
GJ = EI = 1e3
ρA,ρI = 1, 10
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)])

# BCs
τ = 5
θ₀ = -1.5
θ = t -> ifelse.(t.<=τ, θ₀*t/τ, θ₀)
p = t -> 4*tan.(θ(t)/4)
support = create_BC(name="support",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
driver = create_BC(name="driver",beam=beam,node=1,conditions=Dict(:p2A=>t->p(t)))

# Model
rotationDrivenRobotArm = create_Model(name="rotationDrivenRobotArm",beams=[beam],BCs=[support,driver])

# Time variables
tf = 9
Δt = 1e-1

# Create and solve the problem
problem = create_DynamicProblem(model=rotationDrivenRobotArm,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][nElem].nodal.u2[1] for i in 1:length(t)]
u3_tip = [states[i][nElem].nodal.u2[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/rotationDrivenRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=([-L,L],[0,L],[0,L]),save=true,savePath=string(relPath,"/rotationDrivenRobotArm_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,u1_tip/L, lw=lw, label="\$u_1/L\$")
plot!(t,u3_tip/L, lw=lw, label="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/rotationDrivenRobotArm_disp.pdf"))

println("Finished rotationDrivenRobotArm.jl")