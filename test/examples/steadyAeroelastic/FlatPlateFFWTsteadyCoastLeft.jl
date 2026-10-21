using AeroBeams

# Flare angle [rad]
flareAngle = deg2rad(20)

# Root pitch angle
θ = deg2rad(5)

# Airspeed [m/s]
U = 40

# Wing geometry
L = 1
nElem = 20
hingeNode = Int(0.25*nElem) + 1
chord = 0.1
normSparPos = 0.5
thickness = 1e-3

# Sectional stiffness properties
EA, GJ, EIy, EIz = 1e7, 10, 10, 3000

# Sectional inertia properties
ρA = 0.5
ρIy, ρIz = ρA*thickness^2/12, ρA*chord^2/12

# Airfoil and aerodynamic surface
airfoil = deepcopy(FlatPlate)
update_all_airfoil_parameters!(airfoil,U/standard_atmosphere(0).a)
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=false)

# Wing beam
beam = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)],p0=transform_to_WM([0;0;θ],from=:E321),surface=surf,hingedNodes=[hingeNode],hingedNodesDoF=[trues(3)])

# BCs
rootClamp = create_BC(name="clamp",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Hinge axis constraint
localHingeAxis = rotation_tensor_E321([flareAngle, 0, 0]) * [0,1,0]
hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,localHingeAxis=localHingeAxis,masterOnLeft=false)

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Model
FlatPlateFFWTsteadyCoastLeft = create_Model(name="FlatPlateFFWTsteadyCoastLeft",beams=[beam],BCs=[rootClamp],gravityVector=[0,0,-standard_atmosphere(0).g],v_A=t->[0,U,0],hingeAxisConstraints=[hingeAxisConstraint])

# Create and solve problem
problem = create_SteadyProblem(model=FlatPlateFFWTsteadyCoastLeft,systemSolver=systemSolver)
solve!(problem)

# Undeformed positions and coordinates
beam = FlatPlateFFWTsteadyCoastLeft.beams[1]
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
hacState = problem.loadstepHistory.hacStates[end][1]
aeroVars = problem.loadstepHistory.aeroVariables[end]
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
α = [aeroVars[e].unsteadyKin.αₑ for e in 1:nElem]
cn = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
ϕHinge = rad2deg(hacState.ϕ)

# Display results
println("Coast angle = $(round(ϕHinge,sigdigits=4)) deg")

## Plots
# Set paths
relPath = "/test/outputs/figures/FlatPlateFFWTsteadyCoastLeft"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=false,plotDistLoads=true,view=(0,0),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/FlatPlateFFWTsteadyCoastLeft_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan seen from x2-plane
plt1 = plot((r_n1.+u1)/x1_n[end], (r_n3.+u3)/x1_n[end], lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position")
display(plt1)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_disp_x2plane.pdf"))

# Normalized deformed wingspan seen from x3-plane
plt12 = plot((r_n1.+u1)/x1_n[end], (r_n2.+u2)/x1_n[end], lw=lw, label=false, xlabel="Normalized spanwise position", ylabel="Normalized in-plane position", xlims=[0, 1])
display(plt12)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_disp_x3plane.pdf"))

# p1
plt3 = plot(x1_n/x1_n[end], p1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt3)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_p1.pdf"))

# p2
plt2 = plot(x1_n/x1_n[end], p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_p2.pdf"))

# p3
plt2 = plot(x1_n/x1_n[end], p3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_3\$")
display(plt2)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_p3.pdf"))

# Bending moment
plt4 = plot(x1_n/x1_n[end], M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_M2.pdf"))

# Angle of attack
plt5 = plot(x1_e/x1_n[end], α*180/π, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]")
display(plt5)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_aoa.pdf"))

# cn
plt6 = plot(x1_e/x1_n[end], cn, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$c_n\$ ")
display(plt6)
savefig(string(absPath,"/FlatPlateFFWTsteadyCoastLeft_cn.pdf"))

println("Finished FlatPlateFFWTsteadyCoastLeft.jl")