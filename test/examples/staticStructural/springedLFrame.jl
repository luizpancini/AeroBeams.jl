using AeroBeams

# Beam 1
L1 = 2
EIy1 = 1e16
nElem1 = 10
beam1 = create_Beam(name="beam1",length=L1,nElements=nElem1,S=[isotropic_stiffness_matrix(EIy=EIy1)])

# Beam 2
L2 = 1
EIy2 = 1e4
nElem2 = 20
beam2 = create_Beam(name="beam2",length=L2,nElements=nElem2,S=[isotropic_stiffness_matrix(EIy=EIy2)],p0=transform_to_WM([0;π/2;0],from=:E321))

# # Spring
# spring = create_Spring(ku=1e4)
# add_spring_to_beams!(beams=[beam1,beam2],spring=spring,nodesLocalIDs=[div(nElem1,2)+1,nElem2+1],p=[0;π/2;0])

# BCs
F = 1e2
clamp = create_BC(name="clamp",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipForce = create_BC(name="tipForce",beam=beam2,node=nElem2+1,conditions=Dict(:F3b=>F))

# Model
springedLFrame = create_Model(name="springedLFrame",beams=[beam1,beam2],BCs=[clamp,tipForce])

# System solver
NR = create_NewtonRaphson(displayStatus=true)

# Create and solve the problem
problem = create_SteadyProblem(model=springedLFrame,systemSolver=NR)
solve!(problem)

# Get nodal arclength positions and displacement over tip beam
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam1.elements[e].geometry.x1_n1,
                  beam1.elements[e].geometry.x1_n2] for e in 1:nElem1)
u3_b = reduce(vcat,[hist.states[end][e].nodal.u1_b[3],
                    hist.states[end][e].nodal.u2_b[3]] for e in 1:nElem1)

# Tip u3 displacement (in local basis)
println("Tip disp = $(u3_b[end])")

## Plots
# Set paths
relPath = "/test/outputs/figures/springedLFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e2,save=true,savePath=string(relPath,"/springedLFrame_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Plots
plt1 = plot(x1/L2, u3_b/L2, lw=lw, label=false, xlabel="\$x_1/L_2\$", ylabel="\$u_3^{+}/L_2\$")
display(plt1)
savefig(string(absPath,"/springedLFrame_disp.pdf"))

println("Finished springedLFrame.jl")