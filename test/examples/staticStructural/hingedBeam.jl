using AeroBeams

# Option for linear solution
linearSolution = true

# Beam 
L = 2
EIy = 1
nElem = 20
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[midElem+1],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = -1
q = (x1,t) -> ifelse.(x1.>=L/2,q₀,0)
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t=>(x1,t)->[0; 0; q(x1,t)]))
roller = create_BC(name="roller",beam=beam,node=1,conditions=Dict(:u3A=>0))
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
hingedBeam = create_Model(name="hingedBeam",beams=[beam],BCs=[roller,clamp])

# System solver
solver = linearSolution ?
         create_LinearSolver() :
         create_NewtonRaphson()

# Create and solve the problem
problem = create_SteadyProblem(model=hingedBeam,systemSolver=solver)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beam
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

# Linear FEM solution (example 5.2.4 of REDDY - An Introduction to the Finite Element Method - 3rd Ed. - 2006)
u3MidRef = q₀*(L/2)^4/(8*EIy)
p2LeftRef = -q₀*(L/2)^4/(8*L/2*EIy)
p2RightRef = q₀*(L/2)^3/(6*EIy)

# AeroBeams solution
u3Mid = u3[2*midElem+1]
p2Left = p2[2*midElem]
p2Right = p2[2*midElem+1]

# Compute relative errors
ϵ_u3Mid = u3Mid/u3MidRef - 1
ϵ_p2Left = p2Left/p2LeftRef - 1
ϵ_p2Right = p2Right/p2RightRef - 1

# Display relative errors
println("Relative errors: u3Mid = $(round(100ϵ_u3Mid,digits=4)) %, p2Left = $(round(100ϵ_p2Left,digits=4)) %, p2Right = $(round(100ϵ_p2Right,digits=4)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/hingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e0,showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/hingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/hingedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/hingedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt3)
savefig(string(absPath,"/hingedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/hingedBeam_M2.pdf"))

println("Finished hingedBeam.jl")