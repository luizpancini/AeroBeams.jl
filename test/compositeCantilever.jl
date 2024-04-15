using AeroBeams, LinearAlgebra, Plots

# Beam
L = 60.0
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
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
F = 1e5
ω = 20
force = create_BC(name="force",beam=beam,node=nElem+1,types=["F3A"],values=[t->F*sin(ω*t)])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
compositeCantilever = Model(name="compositeCantilever",beams=[beam],BCs=[clamp,force])

# Time variables
tf = 2
Δt = 1e-3

# Create and solve the problem
problem = DynamicProblem(model=compositeCantilever,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][end].u_n2[1] for i in 1:length(t)]
u2_tip = [problem.nodalStatesOverTime[i][end].u_n2[2] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]
p1_tip = [problem.nodalStatesOverTime[i][end].p_n2[1] for i in 1:length(t)]
p2_tip = [problem.nodalStatesOverTime[i][end].p_n2[2] for i in 1:length(t)]
p3_tip = [problem.nodalStatesOverTime[i][end].p_n2[3] for i in 1:length(t)]
F1_root = [problem.nodalStatesOverTime[i][1].F_n1[1] for i in 1:length(t)]
F2_root = [problem.nodalStatesOverTime[i][1].F_n1[2] for i in 1:length(t)]
F3_root = [problem.nodalStatesOverTime[i][1].F_n1[3] for i in 1:length(t)]
M1_root = [problem.nodalStatesOverTime[i][1].M_n1[1] for i in 1:length(t)]
M2_root = [problem.nodalStatesOverTime[i][1].M_n1[2] for i in 1:length(t)]
M3_root = [problem.nodalStatesOverTime[i][1].M_n1[3] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
# Tip u1 
plt1 = plot()
plot!(t,u1_tip, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$u_1\$ [m] ", label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_u1.pdf"))
# Tip u2 
plt2 = plot()
plot!(t,u2_tip, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$u_2\$ [m] ", label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_u2.pdf"))
# Tip u3
plt3 = plot()
plot!(t,u3_tip, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [m] ", label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_u3.pdf"))
# Tip p1 
plt4 = plot()
plot!(t,p1_tip, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$p_1\$", label=false)
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_p1.pdf"))
# Tip p2 
plt5 = plot()
plot!(t,p2_tip, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$p_2\$", label=false)
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_p2.pdf"))
# Tip p3
plt6 = plot()
plot!(t,p3_tip, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Tip \$p_3\$", label=false)
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_p3.pdf"))
# Root F1 
plt7 = plot()
plot!(t,F1_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$F_1^*\$ [N] ", label=false)
display(plt7)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_F1.pdf"))
# Root F2
plt8 = plot()
plot!(t,F2_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$F_2^*\$ [N] ", label=false)
display(plt8)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_F2.pdf"))
# Root F3 
plt9 = plot()
plot!(t,F3_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$F_3^*\$ [N] ", label=false)
display(plt9)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_F3.pdf"))
# Root M1 
plt10 = plot()
plot!(t,M1_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$M_1^*\$ [Nm] ", label=false)
display(plt10)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_M1.pdf"))
# Root M2 
plt11 = plot()
plot!(t,M2_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm] ", label=false)
display(plt11)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_M2.pdf"))
# Root M3 
plt12 = plot()
plot!(t,M3_root, c=:black, linewidth=1, xlabel="\$t\$ [s]", ylabel="Root \$M_3^*\$ [Nm] ", label=false)
display(plt12)
savefig(string(pwd(),"/test/outputs/figures/compositeCantilever_M3.pdf"))

println("Finished compositeCantilever.jl")