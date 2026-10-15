using AeroBeams

# Hinge angle [rad]
θ = deg2rad(90)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[false,true,false]])

# Hinge axis constraint
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=[0;1;0],pHValue=4*tan(θ/4))

# BCs
Fₕ = -1
Fₜ = -1
qᵢ = -1
qₒ = -1
q = (x1,t) -> ifelse.(x1.<=L*(hingeNode-1)/nElem,qᵢ,0) .+ ifelse.(x1.>L*(hingeNode-1)/nElem,qₒ,0)
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t=>(x1,t)->[0, 0, q(x1,t)]))
hingeForce = create_BC(name="hingeForce",beam=beam,node=hingeNode,conditions=Dict(:F3A=>Fₕ))
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:F3A=>Fₜ))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
clampedHingedBeam = create_Model(name="clampedHingedBeam",beams=[beam],BCs=[clamp,hingeForce,tipForce],hingeAxisConstraints=[hingeAxisConstraint])

# System solver
NR = create_NLsolver(method=:newton,linesearch=true)

# Create and solve the problem
problem = create_SteadyProblem(model=clampedHingedBeam,systemSolver=NR)
solve!(problem)

# Get outputs
hist = problem.loadstepHistory
u1 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                  hist.states[end][e].nodal.u2[1]] for e in 1:nElem)
u2 = reduce(vcat,[hist.states[end][e].nodal.u1[2],
                  hist.states[end][e].nodal.u2[2]] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)
p1 = reduce(vcat,[hist.states[end][e].nodal.p1[1],
                  hist.states[end][e].nodal.p2[1]] for e in 1:nElem)
p2 = reduce(vcat,[hist.states[end][e].nodal.p1[2],
                  hist.states[end][e].nodal.p2[2]] for e in 1:nElem)
p3 = reduce(vcat,[hist.states[end][e].nodal.p1[3],
                  hist.states[end][e].nodal.p2[3]] for e in 1:nElem)
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)

# Get nodal arclength positions
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                  beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)

# Compute rotation, displacement and balance moment at the hinge node
ϕHinge = rad2deg(hingeAxisConstraint.state.ϕ)
hingeMoment = hingeAxisConstraint.state.hingeMoment
u3HingeNode = u3[2*hingeNode-1]

# Display results
println("Rotation at hinge node = $(round(ϕHinge,digits=2)) deg")
println("Displacement at hinge node = $(round(u3HingeNode,digits=2)) m")
println("Hinge moment = $(round(hingeMoment,digits=2)) N.m")

## Plots
# Set paths
relPath = "/test/outputs/figures/clampedHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1,showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/clampedHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot((r_n1.+u1)/x1[end], (r_n3.+u3)/x1[end], lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/clampedHingedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
display(plt3)
savefig(string(absPath,"/clampedHingedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt4)
savefig(string(absPath,"/clampedHingedBeam_M2.pdf"))

println("Finished clampedHingedBeam.jl")