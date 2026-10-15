using AeroBeams

# This problem was defined by [Wang and Yu](https://doi.org/10.1063/1.4985091)
# Beam
L = 60
nElem = 10
stiffnessMatrix = [ 2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8;
                    1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5;
                    6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6;
                   -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6;
                   -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7;
                   -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8]
inertiaMatrix = [ 258.053        0.0      0.0       0.0  7.07839  -71.6871;
                      0.0    258.053      0.0  -7.07839      0.0       0.0;
                      0.0        0.0  258.053   71.6871      0.0       0.0;
                      0.0   -7.07839  71.6871     48.59      0.0       0.0;
                  7.07839        0.0      0.0       0.0    2.172       0.0;
                 -71.6871        0.0      0.0       0.0      0.0    46.418]
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
F = 1e5
ω = 20
force = create_BC(name="force",beam=beam,node=nElem+1,conditions=Dict(:F3A=>t->F*sin(ω*t)))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
compositeCantilever = create_Model(name="compositeCantilever",beams=[beam],BCs=[clamp,force])

# Time variables
tf = 2.0
Δt = 1e-3

# Create and solve the problem
problem = create_DynamicProblem(model=compositeCantilever,finalTime=tf,Δt=Δt,displayProgress=false)
@time solve!(problem)

# problem = create_DynamicProblem(model=compositeCantilever,finalTime=tf,Δt=Δt,displayProgress=false)
# @profview_allocs solve!(problem) sample_rate = 0.01

# using Profile
# Profile.clear()
# Profile.init(n=10^7,delay=1e-5)
# problem = create_DynamicProblem(model=compositeCantilever,finalTime=tf,Δt=Δt,displayProgress=false)
# @profview solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][end].nodal.u2[1] for i in 1:length(t)]
u2_tip = [states[i][end].nodal.u2[2] for i in 1:length(t)]
u3_tip = [states[i][end].nodal.u2[3] for i in 1:length(t)]
p1_tip = [states[i][end].nodal.p2[1] for i in 1:length(t)]
p2_tip = [states[i][end].nodal.p2[2] for i in 1:length(t)]
p3_tip = [states[i][end].nodal.p2[3] for i in 1:length(t)]
F1_root = [states[i][1].nodal.F1[1] for i in 1:length(t)]
F2_root = [states[i][1].nodal.F1[2] for i in 1:length(t)]
F3_root = [states[i][1].nodal.F1[3] for i in 1:length(t)]
M1_root = [states[i][1].nodal.M1[1] for i in 1:length(t)]
M2_root = [states[i][1].nodal.M1[2] for i in 1:length(t)]
M3_root = [states[i][1].nodal.M1[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/compositeCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1,plotFrequency=10,plotLimits=([0,L],[-5,5],[-5,5]),save=true,savePath=string(relPath,"/compositeCantilever_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip u1
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_1\$ [m] ")
plot!(t,u1_tip, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/compositeCantilever_u1.pdf"))

# Tip u2 
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_2\$ [m] ")
plot!(t,u2_tip, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/compositeCantilever_u2.pdf"))

# Tip u3
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [m] ")
plot!(t,u3_tip, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/compositeCantilever_u3.pdf"))

# Tip p1 
plt4 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$p_1\$")
plot!(t,p1_tip, c=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/compositeCantilever_p1.pdf"))

# Tip p2 
plt5 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$p_2\$")
plot!(t,p2_tip, c=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/compositeCantilever_p2.pdf"))

# Tip p3
plt6 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$p_3\$")
plot!(t,p3_tip, c=:black, lw=lw, label=false)
display(plt6)
savefig(string(absPath,"/compositeCantilever_p3.pdf"))

# Root F1 
plt7 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_1^*\$ [N]")
plot!(t,F1_root, c=:black, lw=lw, label=false)
display(plt7)
savefig(string(absPath,"/compositeCantilever_F1.pdf"))

# Root F2
plt8 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_2^*\$ [N]")
plot!(t,F2_root, c=:black, lw=lw, label=false)
display(plt8)
savefig(string(absPath,"/compositeCantilever_F2.pdf"))

# Root F3 
plt9 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_3^*\$ [N]")
plot!(t,F3_root, c=:black, lw=lw, label=false)
display(plt9)
savefig(string(absPath,"/compositeCantilever_F3.pdf"))

# Root M1 
plt10 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_1^*\$ [Nm]")
plot!(t,M1_root, c=:black, lw=lw, label=false)
display(plt10)
savefig(string(absPath,"/compositeCantilever_M1.pdf"))

# Root M2 
plt11 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm]")
plot!(t,M2_root, c=:black, lw=lw, label=false)
display(plt11)
savefig(string(absPath,"/compositeCantilever_M2.pdf"))

# Root M3 
plt12 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_3^*\$ [Nm]")
plot!(t,M3_root, c=:black, lw=lw, label=false)
display(plt12)
savefig(string(absPath,"/compositeCantilever_M3.pdf"))

println("Finished compositeCantilever.jl")