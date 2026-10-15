using AeroBeams

# Option for linear solution
linearSolution = true

# Beam
a = b = 1
L = a + b
EIy = 1
nElem = 20
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[midElem+1],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = -1
F₀ = -1
q = (x1,t) -> ifelse.(x1.>=L/2,q₀,0)
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t => (x1,t)->[0; 0; q(x1,t)]))
force = create_BC(name="force",beam=beam,node=midElem+1,conditions=Dict(:F3A=>F₀))
clamp1 = create_BC(name="clamp1",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
biclampedHingedBeam = create_Model(name="biclampedHingedBeam",beams=[beam],BCs=[clamp1,clamp2,force])

# System solver (for nonlinear solution)
solver = linearSolution ?
         create_LinearSolver() :
         create_NewtonRaphson() 

# Create and solve the problem
problem = create_SteadyProblem(model=biclampedHingedBeam,systemSolver=solver)
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

# Linear FEM solution (problem 5.20 of REDDY - An Introduction to the Finite Element Method - 3rd Ed. - 2006)
u3MidRef = 1/EIy*a^3/(a^3+b^3)*(q₀*b^4/8+F₀*b^3/3)
p2RightRef = 1/EIy*b^3/(a^3+b^3)*(q₀*(8*a^3-b^3)/48+F₀*a^3/(2*b))

# The solution for the angle of rotation at the left of the hinge is not given by Reddy, but it is close to the following (verified using a linear beam FE code):
p2LeftRef = -1/EIy*b^3/(a^3+b^3)*(q₀*(8*a^3+b^3)/48+F₀*a^3/(2*b))

# AeroBeams solution
u3Mid = u3[2*midElem+1]
p2Left = p2[2*midElem]
p2Right = p2[2*midElem+1]

# Compute relative errors
ϵ_u3Mid = 1 - u3Mid/u3MidRef
ϵ_p2Left = 1 - p2Left/p2LeftRef
ϵ_p2Right = 1 - p2Right/p2RightRef

# Display relative errors
println("Relative errors: u3Mid = $(round(100*ϵ_u3Mid, sigdigits=2)) %, p2Left = $(round(100*ϵ_p2Left, sigdigits=2)) %, p2Right = $(round(100*ϵ_p2Right, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/biclampedHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
scale = 1
deformationPlot = plot_steady_deformation(problem,scale=scale,showScale=true,scalePos=[0.3,0.1],save=true,savePath=string(relPath,"/biclampedHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/biclampedHingedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/biclampedHingedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt3)
savefig(string(absPath,"/biclampedHingedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/biclampedHingedBeam_M2.pdf"))

println("Finished biclampedHingedBeam.jl")