using AeroBeams

# Beam
L,b,H = 1, 1e-2, 1e-3
E,ρ = 200e9, 7.8e3
A,Iy = b*H, b*H^3/12
∞ = 1e12
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)],I=[inertia_matrix(ρA=ρ*A)])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
cantileverUnderSelfWeight = create_Model(name="cantileverUnderSelfWeight",beams=[beam],BCs=[clamp],gravityVector=[0,0,-9.81])

# Set system solver
σ0 = 0
σstep = 0.05
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=cantileverUnderSelfWeight,systemSolver=NR)
@time solve!(problem)

# problem = create_SteadyProblem(model=cantileverUnderSelfWeight,systemSolver=NR)
# @profview_allocs solve!(problem) sample_rate=0.1

# Get nodal arclength positions, displacements and forces over the beams
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u1 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                  hist.states[end][e].nodal.u2[1]] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)

## Plots
# Set paths
relPath = "/test/outputs/figures/cantileverUnderSelfWeight"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/cantileverUnderSelfWeight_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u1
plt0 = plot(x1/L, u1/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_1/L\$")
display(plt0)
savefig(string(absPath,"/cantileverUnderSelfWeight_u1.pdf"))

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/cantileverUnderSelfWeight_u3.pdf"))

# F3
plt2 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/cantileverUnderSelfWeight_F3.pdf"))

# M2
plt3 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/cantileverUnderSelfWeight_M2.pdf"))

println("Finished cantileverUnderSelfWeight.jl")