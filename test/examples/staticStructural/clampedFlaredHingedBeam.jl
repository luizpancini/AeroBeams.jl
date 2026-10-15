using AeroBeams

# Hinge flare angle [rad]
Λ = deg2rad(30)

# Beam rotation angle about axis a3 [rad]
α = deg2rad(0)

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
surf = create_AeroSurface(airfoil=deepcopy(NACA0018),c=0.1,normSparPos=0.5)

# Gravity
g = standard_atmosphere(0).g

# Beam
L = 1
EIy = 1e1
ρA = 1e-1
nElem = 10
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],surface=surf,p0=transform_to_WM([α;0;0], from=:E321))

# Hinge axis constraint (with hinge axis defined in the local, undeformed beam basis)
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=rotation_tensor_E321([-Λ; 0; 0]) * AeroBeams.a2,pHValue=4*tan(π/2/4))

# Springs around hinge
kT = kOOP = kIP = 1e-6
springT = create_Spring(kp=kT)
springOOP = create_Spring(kp=kOOP)
springIP = create_Spring(kp=kIP)
add_spring_to_beams!(beams=[beam,beam],nodesLocalIDs=[hingeNode-1,hingeNode+1],spring=springT)
add_spring_to_beams!(beams=[beam,beam],nodesLocalIDs=[hingeNode-1,hingeNode+1],spring=springOOP,p=[π/2;0;0])
add_spring_to_beams!(beams=[beam,beam],nodesLocalIDs=[hingeNode-1,hingeNode+1],spring=springIP,p=[0;π/2;0])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# System solver
NR = create_NLsolver(method=:newton,linesearch=true)

# Model
clampedFlaredHingedBeam = create_Model(name="clampedFlaredHingedBeam",beams=[beam],BCs=[clamp],gravityVector=[0,0,-g],hingeAxisConstraints=[hingeAxisConstraint])

# Create and solve the problem
problem = create_SteadyProblem(model=clampedFlaredHingedBeam,systemSolver=NR)
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

# Extract rotation angle and hinge moment about the hinge node
ϕHinge = rad2deg(hingeAxisConstraint.state.ϕ)
hingeMoment = hingeAxisConstraint.state.hingeMoment

# Display results
println("Rotation at hinge node: $(round(ϕHinge,digits=2)) deg")
println("Hinge moment: $(round(hingeMoment,digits=2)) N.m")

## Plots
# Set paths
relPath = "/test/outputs/figures/clampedFlaredHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,15),scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/clampedFlaredHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Position on x1-x2 plane
plt1 = plot((r_n1.+u1)/L, (r_n2.+u2)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_2/L\$")
display(plt1)
savefig(string(absPath,"/clampedFlaredHingedBeam_u2.pdf"))

# Position on x1-x3 plane
plt1 = plot((r_n1.+u1)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedFlaredHingedBeam_u3.pdf"))

# p1
plt2 = plot(x1/L, p1, lw=lw, label=false, xlabel="\$s/L\$", ylabel="\$p_1\$")
display(plt2)
savefig(string(absPath,"/clampedFlaredHingedBeam_p1.pdf"))

# p2
plt3 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$s/L\$", ylabel="\$p_2\$")
display(plt3)
savefig(string(absPath,"/clampedFlaredHingedBeam_p2.pdf"))

# p3
plt4 = plot(x1/L, p3, lw=lw, label=false, xlabel="\$s/L\$", ylabel="\$p_3\$")
display(plt4)
savefig(string(absPath,"/clampedFlaredHingedBeam_p3.pdf"))

# F3
plt5 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$s/L\$", ylabel="\$F_3^{\\star}\$ [N]")
display(plt5)
savefig(string(absPath,"/clampedFlaredHingedBeam_F3.pdf"))

# M2
plt6 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$s/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt6)
savefig(string(absPath,"/clampedFlaredHingedBeam_M2.pdf"))

println("Finished clampedFlaredHingedBeam.jl")