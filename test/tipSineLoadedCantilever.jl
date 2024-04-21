using AeroBeams, LinearAlgebra, Plots

# Beam
L,b,H = 1.0,0.1,0.1
E,ρ = 200e9,7.8e3
A,Iy,Iz, = b*H,b*H^3/12,H*b^3/12
Is = Iy+Iz
∞ = 1e12
nElements = 20
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*Iy,∞])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElements,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
F = 1e3
ω = 200
force = create_BC(name="force",beam=beam,node=nElements+1,types=["F3A"],values=[t->F*sin(ω*t)])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
tipSineLoadedCantilever = create_Model(name="tipSineLoadedCantilever",beams=[beam],BCs=[clamp,force])

# Time variables
T = 2*π/ω
tf = 20*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=tipSineLoadedCantilever,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]
F3_root = [problem.nodalStatesOverTime[i][1].F_n1[3] for i in 1:length(t)]
M2_root = [problem.nodalStatesOverTime[i][1].M_n1[2] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
# Tip u3
plt3 = plot()
plot!(t,u3_tip*1e3, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [mm] ", label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/tipSineLoadedCantilever_u3.pdf"))
# Root F3 
plt9 = plot()
plot!(t,F3_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$F_3^*\$ [N] ", label=false)
display(plt9)
savefig(string(pwd(),"/test/outputs/figures/tipSineLoadedCantilever_F3.pdf"))
# Root M2 
plt11 = plot()
plot!(t,M2_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm] ", label=false)
display(plt11)
savefig(string(pwd(),"/test/outputs/figures/tipSineLoadedCantilever_M2.pdf"))

println("Finished tipSineLoadedCantilever.jl")