using AeroBeams

# Hinge node
hingeNode = 12

# Flare angle [rad]
flareAngle = deg2rad(20)

# Root pitch angle
θ = deg2rad(7)

# Airspeed [m/s]
U = 40

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true,displayStatus=true)

# Model
PazyFFWTsteadyCoast,nElem,L,_ = create_Pazy(hinged=true,hingeNode=hingeNode,flareAngle=flareAngle,airspeed=U,pitchAngle=θ)

# Create and solve problem
problem = create_SteadyProblem(model=PazyFFWTsteadyCoast,systemSolver=systemSolver)
solve!(problem)

# Undeformed positions and coordinates
beam = PazyFFWTsteadyCoast.beams[1]
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
relPath = "/test/outputs/figures/PazyFFWTsteadyCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=false,plotDistLoads=true,view=(0,0),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteadyCoast_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan seen from x2-plane
plt1 = plot((r_n1.+u1)/x1_n[end], (r_n3.+u3)/x1_n[end], lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position")
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyCoast_disp_x2plane.pdf"))

# Normalized deformed wingspan seen from x3-plane
plt12 = plot((r_n1.+u1)/x1_n[end], (r_n2.+u2)/x1_n[end], lw=lw, label=false, xlabel="Normalized spanwise position", ylabel="Normalized in-plane position", xlims=[0, 1])
display(plt12)
savefig(string(absPath,"/PazyFFWTsteadyCoast_disp_x3plane.pdf"))

# p1
plt3 = plot(x1_n/x1_n[end], p1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt3)
savefig(string(absPath,"/PazyFFWTsteadyCoast_p1.pdf"))

# p2
plt2 = plot(x1_n/x1_n[end], p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyCoast_p2.pdf"))

# p3
plt2 = plot(x1_n/x1_n[end], p3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_3\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyCoast_p3.pdf"))

# Bending moment
plt4 = plot(x1_n/x1_n[end], M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/PazyFFWTsteadyCoast_M2.pdf"))

# Angle of attack
plt5 = plot(x1_e/x1_n[end], α*180/π, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]")
display(plt5)
savefig(string(absPath,"/PazyFFWTsteadyCoast_aoa.pdf"))

# cn
plt6 = plot(x1_e/x1_n[end], cn, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$c_n\$ ")
display(plt6)
savefig(string(absPath,"/PazyFFWTsteadyCoast_cn.pdf"))

println("Finished PazyFFWTsteadyCoast.jl")