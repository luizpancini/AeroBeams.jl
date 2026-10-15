using AeroBeams

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 12
foldAngle = -deg2rad(10)
flareAngle = deg2rad(20)

# Root pitch angle
θ = deg2rad(7)

# Airspeed [m/s]
U = 40

# Pazy wing with flared folding tip
PazyFFWTsteadyFixedFold,nElem,L,_ = create_Pazy(hinged=true,hingeNode=hingeNode,foldAngle=foldAngle,flareAngle=flareAngle,airspeed=U,pitchAngle=θ)

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true,displayStatus=true)

# Create and solve problem
problem = create_SteadyProblem(model=PazyFFWTsteadyFixedFold,systemSolver=systemSolver)
solve!(problem)

# Undeformed positions and coordinates
beam = PazyFFWTsteadyFixedFold.beams[1]
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
hingeMoment = hacState.hingeMoment
u3HingeNode = u3[2*hingeNode-1]
u3Tip = u3[end]

# Display results
println("Total rotation at hinge node = $(round(ϕHinge,sigdigits=4)) deg")
println("Moment necessary to impose the hinge angle = $(round(hingeMoment,sigdigits=4)) N.m")
println("OOP displacements: at hinge node = $(round(u3HingeNode,sigdigits=4)) m, at tip = $(round(u3Tip,sigdigits=4)) m")

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyFixedFold"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,view=(60,30),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteadyFixedFold_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# OOP displacement
plt1 = plot((r_n1.+u1)/x1_n[end], (r_n3.+u3)/x1_n[end], lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position")
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_u3.pdf"))

# p1
plt3 = plot(x1_n/x1_n[end], p1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt3)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_p1.pdf"))

# p2
plt2 = plot(x1_n/x1_n[end], p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_p2.pdf"))

# p3
plt2 = plot(x1_n/x1_n[end], p3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_3\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_p3.pdf"))

# Bending moment
plt4 = plot(x1_n/x1_n[end], M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_M2.pdf"))

# Angle of attack
plt5 = plot(x1_e/x1_n[end], α*180/π, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]")
display(plt5)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_aoa.pdf"))

# cn
plt6 = plot(x1_e/x1_n[end], cn, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$c_n\$ ")
display(plt6)
savefig(string(absPath,"/PazyFFWTsteadyFixedFold_cn.pdf"))

println("Finished PazyFFWTsteadyFixedFold.jl")