using AeroBeams

# Beam rotation angle about a3-axis
α = deg2rad(30)

# Beam rotation angle about a2'-axis
β = deg2rad(-0)

# Hinge angle [rad]
θ = deg2rad(-45)

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
surf = create_AeroSurface(airfoil=deepcopy(NACA0018),c=0.1,normSparPos=0.25)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],p0=transform_to_WM([α;β;0];from=:E321),surface=surf)

# Hinge axis constraint
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=[0;1;0],pHValue=4*tan(θ/4))

# BCs
Fₕ = -1
Fₜ = -1
qᵢ = -1
qₒ = -1
q = (x1,t) -> ifelse.(x1.<=L*(hingeNode-1)/nElem,qᵢ,0) .+ ifelse.(x1.>L*(hingeNode-1)/nElem,qₒ,0)
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t=>(x1,t)->[0; 0; q(x1,t)]))
hingeForce = create_BC(name="hingeForce",beam=beam,node=hingeNode,conditions=Dict(:F3A=>Fₕ))
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:F3A=>Fₜ))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
clampedHingedBeamRotated = create_Model(name="clampedHingedBeamRotated",beams=[beam],BCs=[clamp,hingeForce,tipForce],hingeAxisConstraints=[hingeAxisConstraint])

# System solver
NR = create_NLsolver(method=:newton,linesearch=true)

# Create and solve the problem
problem = create_SteadyProblem(model=clampedHingedBeamRotated,systemSolver=NR)
solve!(problem)

# Get outputs
hist = problem.loadstepHistory
u1 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                  hist.states[end][e].nodal.u2[1]] for e in 1:nElem)
u2 = reduce(vcat,[hist.states[end][e].nodal.u1[2],
                  hist.states[end][e].nodal.u2[2]] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)
p1_b = reduce(vcat,[hist.states[end][e].nodal.p1_b[1],
                  hist.states[end][e].nodal.p2_b[1]] for e in 1:nElem)
p2_b = reduce(vcat,[hist.states[end][e].nodal.p1_b[2],
                  hist.states[end][e].nodal.p2_b[2]] for e in 1:nElem)
p3_b = reduce(vcat,[hist.states[end][e].nodal.p1_b[3],
                  hist.states[end][e].nodal.p2_b[3]] for e in 1:nElem)
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
u3Tip = u3[end]

# Display results
println("Rotation at hinge node = $(round(ϕHinge,digits=2)) deg")
println("Displacement at hinge node = $(round(u3HingeNode,digits=2)) m")
println("Displacement at tip = $(round(u3Tip,digits=2)) m")
println("Hinge moment = $(round(hingeMoment,digits=2)) N.m")

## Plots
# Set paths
relPath = "/test/outputs/figures/clampedHingedBeamRotated"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,15),showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/clampedHingedBeamRotated_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3 (seen from x2-plane)
plt1 = plot((r_n1.+u1)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedBeamRotated_u3_x2plane.pdf"))

# u3 (seen from x1-plane)
plt1 = plot((r_n2.+u2)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_2/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedBeamRotated_u3_x1plane.pdf"))

# p2_b (local basis)
plt2 = plot(x1/L, p2_b, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$p_2^+\$")
display(plt2)
savefig(string(absPath,"/clampedHingedBeamRotated_p2_b.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$F_3^{\\star}\$ [N]")
display(plt3)
savefig(string(absPath,"/clampedHingedBeamRotated_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt4)
savefig(string(absPath,"/clampedHingedBeamRotated_M2.pdf"))

println("Finished clampedHingedBeamRotated.jl")