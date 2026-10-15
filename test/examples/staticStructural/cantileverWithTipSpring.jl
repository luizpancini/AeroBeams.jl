using AeroBeams

# Beam
L = 1
EIy = 1e4
θ = deg2rad(30)
nElem = 10
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],p0=transform_to_WM([0;θ;0]; from=:E321))

# Set spring in perpendicular bending direction at tip node (rotation of π/2 around x2-axis), with a custom force function and a pre-deformation
f = (k,x) -> k * x.^1 # + 100 * k * x.^3
Δu0 = 0.0
spring = create_Spring(ku=1e4,forceFunction=f)
add_spring_to_beam!(beam=beam,spring=spring,nodeID=nElem+1,p=[0;π/2;0],Δu0=Δu0)

# BCs
F = 1e3
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:F3b=>F))

# Model
cantileverWithTipSpring = create_Model(name="cantileverWithTipSpring",beams=[beam],BCs=[clamp,tipForce])

# Create and solve the problem
problem = create_SteadyProblem(model=cantileverWithTipSpring)
@time solve!(problem)
# @profview_allocs solve!(problem) sample_rate=1

# Get nodal arclength positions, displacement and loads over the beam
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1_b[3],
                  hist.states[end][e].nodal.u2_b[3]] for e in 1:nElem)
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)

# Tip displacement
println("Tip disp = $(u3[end])")

## Plots
# Set paths
relPath = "/test/outputs/figures/cantileverWithTipSpring"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/cantileverWithTipSpring_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/cantileverWithTipSpring_u3.pdf"))

# F3
plt2 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/cantileverWithTipSpring_F3.pdf"))

# M2
plt3 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/cantileverWithTipSpring_M2.pdf"))

println("Finished cantileverWithTipSpring.jl")