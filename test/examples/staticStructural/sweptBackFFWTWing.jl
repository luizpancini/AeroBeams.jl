using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(-30)

# Flare angle [rad]
Γ = deg2rad(20)

# Fold angle [rad]
θ = deg2rad(-90)

# Hinge axis (defined in the local, undeformed beam basis)
localHingeAxis = rotation_tensor_E321([-Γ, 0, 0]) * AeroBeams.a2

# Aerodynamic surface (for better visualization only, no aerodynamics involved)
airfoil = deepcopy(NACA0018)
chord = 0.1
normSparPos = 0.25
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = 3*div(nElem,4)+1
outboardElem = hingeNode
inboardElem = outboardElem-1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[true,true,true]],p0=transform_to_WM([Λ,0,0],from=:E321),surface=surf)

# Hinge axis constraint
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=localHingeAxis,pHValue=4*tan(θ/4))

# BCs
Fₜ = 1e-1
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,conditions=Dict(:Ff3A=>Fₜ))

# Model
sweptBackFFWTWing = create_Model(name="sweptBackFFWTWing",beams=[beam],BCs=[clamp,tipFollowerForce],hingeAxisConstraints=[hingeAxisConstraint])

# System solver
systemSolver = create_NLsolver(displayStatus=true,method=:newton,linesearch=true)

# Create and solve the problem
problem = create_SteadyProblem(model=sweptBackFFWTWing,systemSolver=systemSolver)
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

# Fold angle and hinge moment
c = hingeAxisConstraint.state
ϕHinge = rad2deg(c.ϕ)
hingeMoment = c.hingeMoment

# Display results
println("Hinge rotation angle = $(round(ϕHinge,digits=1)) deg")
println("Hinge moment = $(round(hingeMoment,sigdigits=4))")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptBackFFWTWing"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,view=(15,30),save=true,savePath=string(relPath,"/sweptBackFFWTWing_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Bending moment
plt_M2 = plot(x1_n,M2,xlabel="\$x_1\$",ylabel="\$M_2\$ [N.m]",lw=lw,label=false)
display(plt_M2)
savefig(string(absPath,"/sweptBackFFWTWing_M2.pdf"))

println("Finished sweptBackFFWTWing.jl")