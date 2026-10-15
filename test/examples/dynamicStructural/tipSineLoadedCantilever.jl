using AeroBeams

# Beam
L,b,H = 1.0, 0.1, 0.1
E,ρ = 200e9, 7.8e3
A,Iy,Iz, = b*H, b*H^3/12, H*b^3/12
∞ = 1e12
nElements = 20
beam = create_Beam(name="beam",length=L,nElements=nElements,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)])

# BCs
F = 1e3
ω = 200
force = create_BC(name="force",beam=beam,node=nElements+1,conditions=Dict(:F3A=>t->F*sin(ω*t)))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
tipSineLoadedCantilever = create_Model(name="tipSineLoadedCantilever",beams=[beam],BCs=[clamp,force])

# Time variables
T = 2π/ω
cycles = 20
tf = cycles * T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=tipSineLoadedCantilever,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u3_tip = [states[i][end].nodal.u2[3] for i in 1:length(t)]
F3_root = [states[i][1].nodal.F1[3] for i in 1:length(t)]
M2_root = [states[i][1].nodal.M1[2] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/tipSineLoadedCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1e3,plotFrequency=10,fps=20,plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),save=true,savePath=string(relPath,"/tipSineLoadedCantilever_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
gr()
lw = 1

# Tip u3
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [mm]")
plot!(t,u3_tip*1e3, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/tipSineLoadedCantilever_u3.pdf"))

# Root F3 
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_3^*\$ [N]")
plot!(t,F3_root, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/tipSineLoadedCantilever_F3.pdf"))

# Root M2 
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm]")
plot!(t,M2_root, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/tipSineLoadedCantilever_M2.pdf"))

println("Finished tipSineLoadedCantilever.jl")