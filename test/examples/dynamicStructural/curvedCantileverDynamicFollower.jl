using AeroBeams

# Beam
R,θ = 100, π/4
L = R*θ
EA,GA,GJ,EI = 1e7, 5e6, 833_333, 833_333
ρA,ρI = 1, 10
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],k=[0; 0; 1/R])

# BCs
F₀ = 100
F = t -> F₀*t
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:Ff3A=>(t)->F(t)))

# Model
curvedCantileverDynamicFollower = create_Model(name="curvedCantileverDynamicFollower",beams=[beam],BCs=[clamp,tipForce])

# Time variables
tf = 10
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=curvedCantileverDynamicFollower,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][end].nodal.u2[1] for i in 1:length(t)]
u2_tip = [states[i][end].nodal.u2[2] for i in 1:length(t)]
u3_tip = [states[i][end].nodal.u2[3] for i in 1:length(t)]
F1_root = [states[i][1].nodal.F1[1] for i in 1:length(t)]
F2_root = [states[i][1].nodal.F1[2] for i in 1:length(t)]
F3_root = [states[i][1].nodal.F1[3] for i in 1:length(t)]
M1_root = [states[i][1].nodal.M1[1] for i in 1:length(t)]
M2_root = [states[i][1].nodal.M1[2] for i in 1:length(t)]
M3_root = [states[i][1].nodal.M1[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/curvedCantileverDynamicFollower"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1,plotFrequency=4,plotLimits=([-L,L],[-L,L],[-L,L]),save=true,savePath=string(relPath,"/curvedCantileverDynamicFollower_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
dispLabels=["\$u_1/R\$" "\$u_2/R\$" "\$u_3/R\$"]
forceLabels=["\$F_1^*/F_0\$" "\$F_2^*/F_0\$" "\$F_3^*/F_0\$"]
momLabels=["\$M_1^*/(F_0R)\$" "\$M_2^*/(F_0R)\$" "\$M_3^*/(F_0R)\$"]
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip displacements")
plot!(t,[u1_tip/R,u2_tip/R,u3_tip/R], lw=lw, label=dispLabels)
display(plt1)
savefig(string(absPath,"/curvedCantileverDynamicFollower_disp.pdf"))

# Nomalized root forces
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Root forces")
plot!(t,[F1_root/F₀,F2_root/F₀,F3_root/F₀], lw=lw, label=forceLabels)
display(plt2)
savefig(string(absPath,"/curvedCantileverDynamicFollower_forces.pdf"))

# Nomalized root moments
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root moments")
plot!(t,[M1_root/(F₀*R),M2_root/(F₀*R),M3_root/(F₀*R)], lw=lw, label=momLabels)
display(plt3)
savefig(string(absPath,"/curvedCantileverDynamicFollower_moments.pdf"))

println("Finished curvedCantileverDynamicFollower.jl")