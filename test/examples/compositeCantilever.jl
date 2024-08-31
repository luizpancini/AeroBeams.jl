using AeroBeams, LinearAlgebra

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
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
F = 1e5
ω = 20
force = create_BC(name="force",beam=beam,node=nElem+1,types=["F3A"],values=[t->F*sin(ω*t)])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
compositeCantilever = create_Model(name="compositeCantilever",beams=[beam],BCs=[clamp,force])

# Time variables
tf = 2
Δt = 1e-3

# Create and solve the problem
problem = create_DynamicProblem(model=compositeCantilever,finalTime=tf,Δt=Δt)
solve!(problem)

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

println("Finished compositeCantilever.jl")