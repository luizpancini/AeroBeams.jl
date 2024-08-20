using AeroBeams, LinearAlgebra, Plots

# Gravity
g = 9.81

# Initial angles of release
θ₀ = π/2

# Beam 
L,b,H = 2,0.01,0.01
A,I = b*H,b*H^3/12
E,ρ = 200e9,7.8e3
∞ = 1e16
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*I,∞])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,2*ρ*I,ρ*I,ρ*I])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,(π/2-θ₀),0],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

# BCs
support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])

# Model
doublePendulum = create_Model(name="doublePendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g])

# Time variables
tf = 1
Δt = tf/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=doublePendulum,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[1] for i in 1:length(t)]
u3_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[3] for i in 1:length(t)]
u1_tip = [problem.nodalStatesOverTime[i][end].u_n2[1] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
lw = 2
labels = ["Hinge" "Tip"]
colors = [:black,:blue]
relPath = "/test/outputs/figures/doublePendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,plotFrequency=10,fps=60,scale=1,plotUndeformed=false,plotLimits=[(-L,L),(-L,0),(-L,L)],save=true,savePath=string(relPath,"/doublePendulum_deformation.gif"),displayProgress=true)
# Normalized tip u1 displacement
gr()
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