using AeroBeams

# Gravity
g = standard_atmosphere(0).g

# Initial angle of release
θ₀ = deg2rad(90)

# Beam
L,b,H = 2, 0.01, 0.01
A,I = b*H, b*H^3/12
E,ρ = 200e9, 7.8e3
∞ = 1e16
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*I,∞=∞)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*I,ρIz=ρ*I)],p0=transform_to_WM([0,(π/2-θ₀),0],from=:E321),hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

# BCs
support = create_BC(name="support",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))

# Model
doublePendulum = create_Model(name="doublePendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g])

# Time variables
tf = 1
Δt = tf/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=doublePendulum,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_hinge = [states[i][div(nElem,2)].nodal.u2[1] for i in 1:length(t)]
u3_hinge = [states[i][div(nElem,2)].nodal.u2[3] for i in 1:length(t)]
u1_tip = [states[i][end].nodal.u2[1] for i in 1:length(t)]
u3_tip = [states[i][end].nodal.u2[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/doublePendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=10,fps=60,scale=1,plotUndeformed=false,plotLimits=([-L,L],[-L,0],[-L,L]),save=true,savePath=string(relPath,"/doublePendulum_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
labels = ["Hinge" "Tip"]
colors = [:black,:blue]
gr()

# Normalized tip u1 displacement
plt1 = plot(palette=colors, xlabel="\$t\$ [s]", ylabel="\$u_1/L\$ ")
plot!(t,[u1_hinge,u1_tip]/L, lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/doublePendulum_u1.pdf"))

# Normalized tip u3 displacement
plt2 = plot(palette=colors, xlabel="\$t\$ [s]", ylabel="\$u_3/L\$ ")
plot!(t,[u3_hinge,u3_tip]/L, lw=lw, label=labels)
display(plt2)
savefig(string(absPath,"/doublePendulum_u3.pdf"))

println("Finished doublePendulum.jl")