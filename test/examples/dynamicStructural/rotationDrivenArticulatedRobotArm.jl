using AeroBeams

# Beam
L = 10
EA = GA = 1e6
GJ = EI = 1e5
ρA = ρI = 1
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

# BCs 
τ = 0.5
θ₀ = -1.5
p₀ = 4*tan(θ₀/4)
p = t -> ifelse.(t.<=τ, p₀*t/τ, p₀)
support = create_BC(name="support",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
driver = create_BC(name="driver",beam=beam,node=1,conditions=Dict(:p2A=>t->p(t)))

# Model
rotationDrivenArticulatedRobotArm = create_Model(name="rotationDrivenArticulatedRobotArm",beams=[beam],BCs=[support,driver])

# Time variables
tf = 2.2
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=rotationDrivenArticulatedRobotArm,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][nElem].nodal.u2[1] for i in 1:length(t)]
u3_tip = [states[i][nElem].nodal.u2[3] for i in 1:length(t)]
u1_hinge = [states[i][div(nElem,2)].nodal.u2[1] for i in 1:length(t)]
u3_hinge = [states[i][div(nElem,2)].nodal.u2[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/rotationDrivenArticulatedRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=([-L/2,L],[-L/2,L/2],[-L/10,L]),save=true,savePath=string(relPath,"/rotationDrivenArticulatedRobotArm_deformation.gif"),displayProgress=true)
display(anim)

# Snapshots
plt_snap = plot_snapshots(problem,plotBCs=false,plotAxes=false,plotGrid=false,snapshots=vcat(0,0.6:0.1:2.2),plotLimits=([-L/2,L],[-L/2,L/2],[-L/10,L]),save=true,savePath=string(relPath,"/rotationDrivenArticulatedRobotArm_snapshots.pdf"))
display(plt_snap)

# Plot configurations
lw = 2
labels = ["Tip" "Hinge"]
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="\$u_1/L\$")
plot!(t,[u1_tip/L, u1_hinge/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/rotationDrivenArticulatedRobotArm_u1.pdf"))
plt2 = plot(xlabel="\$t\$ [s]", ylabel="\$u_3/L\$")
plot!(t,[u3_tip/L, u3_hinge/L], lw=lw, label=labels)
display(plt2)
savefig(string(absPath,"/rotationDrivenArticulatedRobotArm_u3.pdf"))

println("Finished rotationDrivenArticulatedRobotArm.jl")