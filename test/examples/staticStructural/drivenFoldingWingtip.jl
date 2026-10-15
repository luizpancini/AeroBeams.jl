using AeroBeams

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.25
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam
L = 1
EIy = 1
nElem = 20
midElem = div(nElem,2)
hingedNode = midElem+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingedNode],hingedNodesDoF=[[true,true,true]],surface=surf)

# Hinge fold angle [rad]
foldAngle = deg2rad(135)

# Hinge flare angle and local axis
flareAngle = deg2rad(20)
localHingeAxis = rotation_tensor_E321([-flareAngle; 0; 0]) * AeroBeams.a2

# Hinge axis
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=localHingeAxis,pHValue=4*tan(foldAngle/4))

# BCs
qᵢ = -1e0
q₀ = -1e-1
q = (x1,t) -> ifelse.(x1.<L/2,qᵢ,0) .+ ifelse.(x1.>=L/2,q₀,0)
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t=>(x1,t)->[0, 0, q(x1,t)]))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
drivenFoldingWingtip = create_Model(name="drivenFoldingWingtip",beams=[beam],BCs=[clamp],hingeAxisConstraints=[hingeAxisConstraint])

# System solver
systemSolver = create_NLsolver(displayStatus=true,method=:newton,linesearch=true)

# Create and solve problem
problem = create_SteadyProblem(model=drivenFoldingWingtip,systemSolver=systemSolver)
solve!(problem)

# Undeformed positions and coordinates
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)
x1_n = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                    beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)
x1_e = [beam.elements[e].geometry.x1 for e in 1:nElem]

# Solution outputs
state = problem.loadstepHistory.states[end]
u1 = reduce(vcat,[state[e].nodal.u1[1],
                  state[e].nodal.u2[1]] for e in 1:nElem)
u2 = reduce(vcat,[state[e].nodal.u1[2],
                  state[e].nodal.u2[2]] for e in 1:nElem)                  
u3 = reduce(vcat,[state[e].nodal.u1[3],
                  state[e].nodal.u2[3]] for e in 1:nElem)
p1 = reduce(vcat,[state[e].nodal.p1[1],
                  state[e].nodal.p2[1]] for e in 1:nElem)
p2 = reduce(vcat,[state[e].nodal.p1[2],
                  state[e].nodal.p2[2]] for e in 1:nElem)                  
p3 = reduce(vcat,[state[e].nodal.p1[3],
                  state[e].nodal.p2[3]] for e in 1:nElem)
F3 = reduce(vcat,[state[e].nodal.F1[3],
                  state[e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[state[e].nodal.M1[2],
                  state[e].nodal.M2[2]] for e in 1:nElem)

# Get rotation at the hinge node
c = hingeAxisConstraint.state
pHinge = c.pH
ϕHinge = -rad2deg(c.ϕ)

# Display results
println("Hinge fold angle = $(round(ϕHinge,digits=1)) deg")

## Plots
# Set paths
relPath = "/test/outputs/figures/drivenFoldingWingtip"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Displacements
plt_u1 = plot(x1_n,u1,ylabel="\$u_1\$", lw=lw, label=false)
plt_u2 = plot(x1_n,u2,ylabel="\$u_2\$", lw=lw, label=false)
plt_u3 = plot(x1_n,u3,xlabel="\$x_1\$",ylabel="\$u_3\$", lw=lw, label=false)
plt_u = plot(plt_u1, plt_u2, plt_u3, layout=(3,1))
display(plt_u)
savefig(string(absPath,"/drivenFoldingWingtip_u.pdf"))

# Rotations
plt_p1 = plot(x1_n,p1,ylabel="\$p_1\$", lw=lw, label=false)
plt_p2 = plot(x1_n,p2,ylabel="\$p_2\$", lw=lw, label=false)
plt_p3 = plot(x1_n,p3,xlabel="\$x_1\$",ylabel="\$p_3\$", lw=lw, label=false)
plt_p = plot(plt_p1, plt_p2, plt_p3, layout=(3,1))
display(plt_p)
savefig(string(absPath,"/drivenFoldingWingtip_p.pdf"))

# Forces and moments
plt_F3 = plot(x1_n,F3,ylabel="\$F_3\$ [N]", lw=lw, label=false)
plt_M2 = plot(x1_n,M2,xlabel="\$x_1\$",ylabel="\$M_2\$ [N.m]", lw=lw, label=false)
plt_FM = plot(plt_F3, plt_M2, layout=(2,1))
display(plt_FM)
savefig(string(absPath,"/drivenFoldingWingtip_FM.pdf"))

# 3D position
plt_pos = plot(r_n1+u1, r_n2+u2, r_n3+u3, xlabel="\$x_1\$", ylabel="\$x_2\$", zlabel="\$x_3\$", lw=lw, label=false, aspect_ratio=:equal)
display(plt_pos)
savefig(string(absPath,"/drivenFoldingWingtip_pos.pdf"))

# 2D position (x1-x3 plane)
plt_pos_x1x3 = plot(r_n1+u1, r_n3+u3, xlabel="\$x_1\$", ylabel="\$x_3\$", lw=lw, label=false, aspect_ratio=:equal)
display(plt_pos_x1x3)
savefig(string(absPath,"/drivenFoldingWingtip_pos_x1x3.pdf"))

println("Finished drivenFoldingWingtip.jl")