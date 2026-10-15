using AeroBeams

# Choose to calculate linear or nonlinear solution (they differ more and more as the distributed load q increases)
linearSolution = false

# Beam frame
L = 1
EI = 200e9 * 5e-6
∞ = 1e12
stiffnessMatrix = isotropic_stiffness_matrix(EIy=EI,EIz=EI,∞=∞)
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[stiffnessMatrix],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[true,false,true]])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[stiffnessMatrix],p0=transform_to_WM([π/2;0;0],from=:E321),hingedNodes=[1],hingedNodesDoF=[[true,false,true]])

# BCs
q = 1e5
add_distributed_loads_to_beam!(beam=beam2,loads=Dict(:f_A_of_x1t=>(x1,t)->[q; 0; q]))
pin1 = create_BC(name="pin1",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
pin2 = create_BC(name="pin2",beam=beam1,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
pin3 = create_BC(name="pin3",beam=beam2,node=nElem+1,conditions=Dict(:u1A=>0,:u3A=>0,:p2A=>0))

# Model
hingedTFrame = create_Model(name="hingedTFrame",beams=[beam1,beam2],BCs=[pin1,pin2,pin3],connections=[(beam2=>beam1,1=>div(nElem,2)+1)])

# System solver
solver = linearSolution ?
         create_LinearSolver() :
         create_NewtonRaphson(displayStatus=true)

# Create and solve the problem
problem = create_SteadyProblem(model=hingedTFrame,systemSolver=solver)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
hist = problem.loadstepHistory
x1_beam1 = reduce(vcat,[beam1.elements[e].geometry.x1_n1,
                        beam1.elements[e].geometry.x1_n2] for e in 1:nElem)
x1_beam2 = reduce(vcat,[beam2.elements[e].geometry.x1_n1,
                        beam2.elements[e].geometry.x1_n2] for e in 1:nElem)
u1_beam1 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                        hist.states[end][e].nodal.u2[1]] for e in 1:nElem)     
u3_beam1 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                        hist.states[end][e].nodal.u2[3]] for e in 1:nElem)
F3_beam1 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                        hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2_beam1 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                        hist.states[end][e].nodal.M2[2]] for e in 1:nElem)
u1_beam2 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                        hist.states[end][e].nodal.u2[1]] for e in nElem+1:2*nElem)
u3_beam2 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                        hist.states[end][e].nodal.u2[3]] for e in nElem+1:2*nElem)

## Plots
# Set paths
relPath = "/test/outputs/figures/hingedTFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=300,view=(60,20),legendPos=(0.4,0.8),save=true,savePath=string(relPath,"/hingedTFrame_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
beamLabels = ["Beam 1" "Beam 2"]
forceLabels = ["\$F_3\$" "\$M_2\$"]
gr()

# u1
plt1 = plot([x1_beam1 x1_beam2], [u1_beam1 u1_beam2]/L, lw=lw, label=beamLabels, xlabel="\$x_1/L\$", ylabel="\$u_1/L\$")
display(plt1)
savefig(string(absPath,"/hingedTFrame_u1.pdf"))

# u3
plt2 = plot([x1_beam1 x1_beam2], [u3_beam1 u3_beam2]/L, lw=lw, label=beamLabels, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt2)
savefig(string(absPath,"/hingedTFrame_u3.pdf"))

# Internal loads
plt3 = plot(x1_beam1, [F3_beam1 M2_beam1], lw=lw, label=forceLabels, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N], \$M_2\$ [N.m]", title="Internal loads on beam 1")
display(plt3)
savefig(string(absPath,"/hingedTFrame_loads.pdf"))

println("Finished hingedTFrame.jl")