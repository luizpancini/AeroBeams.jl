using AeroBeams

# Flag for linear solution
linearSolution = false

# Beam
L = 1.0
EIy = 1e7
∞ = 1e12
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy,∞=∞)])

# BCs
q₀ = 1e6
q = (x1,t) -> q₀*x1
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t=>(x1,t)->[0; 0; q(x1,t)]))
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
roller = create_BC(name="roller",beam=beam,node=1,conditions=Dict(:u3A=>0))

# Model
triangleLoadBeam = create_Model(name="triangleLoadBeam",beams=[beam],BCs=[clamp,roller])

# System solver
solver = linearSolution ? 
         create_LinearSolver() : 
         create_NewtonRaphson(displayStatus=false)

# Create and solve the problem
problem = create_SteadyProblem(model=triangleLoadBeam,systemSolver=solver)
solve!(problem)

# println("Solving with allocation profiler...")
# problem = create_SteadyProblem(model=triangleLoadBeam,systemSolver=solver)
# @profview_allocs solve!(problem) sample_rate=0.1

# println("Solving with time profiler...")
# problem = create_SteadyProblem(model=triangleLoadBeam,systemSolver=solver)
# using Profile
# Profile.clear()
# Profile.init(n=10^7, delay=1e-4)
# @profview solve!(problem)

println("Solving with time measurement:")
problem = create_SteadyProblem(model=triangleLoadBeam,systemSolver=solver)
@time solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)

## Plots
# Set paths
relPath = "/test/outputs/figures/triangleLoadBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e3,save=true,savePath=string(relPath,"/triangleLoadBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Displacement
plt1 = plot(x1/L, u3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3\$ [m]")
display(plt1)
savefig(string(absPath,"/triangleLoadBeam_u3.pdf"))

# Force
plt2 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/triangleLoadBeam_F3.pdf"))

# Moment
plt3 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/triangleLoadBeam_M2.pdf"))

println("Finished triangleLoadBeam.jl")