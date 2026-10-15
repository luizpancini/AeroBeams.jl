using AeroBeams

# Option to solve linear problem
linearSolution = true

# Beam 
L = 2
EA = 1e6
EIy = 1e6
nElem = 20
hingeNode = div(nElem,2) + 1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[false,true,false]])

# Spring
spring = create_Spring(kp=1e5)
add_spring_to_beams!(beams=[beam,beam],spring=spring,nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2;0;0])

# BCs
q₀ = 1e4
q = (x1,t) -> ifelse.(x1.>=L/2,-q₀,0)
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t=>(x1,t)->[0; 0; q(x1,t)]))
pin = create_BC(name="pin",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
hingedSpringedBeam = create_Model(name="hingedSpringedBeam",beams=[beam],BCs=[pin,clamp])

# System solver
solver = linearSolution ?
         create_LinearSolver() :
         create_NewtonRaphson(displayStatus=true)

# Create and solve the problem
problem = create_SteadyProblem(model=hingedSpringedBeam,systemSolver=solver)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)
p2 = reduce(vcat,[hist.states[end][e].nodal.p1[2],
                  hist.states[end][e].nodal.p2[2]] for e in 1:nElem)
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)

# Display midpoint displacement
println("Midpoint u3 = $(u3[hingeNode])")

## Plots
# Set paths
relPath = "/test/outputs/figures/hingedSpringedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e2,showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/hingedSpringedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/hingedSpringedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/hingedSpringedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3/1e3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [kN]")
display(plt3)
savefig(string(absPath,"/hingedSpringedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2/1e3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [kN.m]")
display(plt4)
savefig(string(absPath,"/hingedSpringedBeam_M2.pdf"))

println("Finished hingedSpringedBeam.jl")