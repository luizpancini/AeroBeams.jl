using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Hinge node, hinge angle [rad] and flare angle [deg]
hingeNode = 14
hingeAngle = -1*π/2
flareAngle = 10

# Spring stiffness
kSpring = 0e-1

# Root pitch angle
θ = 7*pi/180

# Airspeed
U = 50

# Gravity
g = 9.80665

# Pazy wing with flared folding tip
pazyFFWT,_ = create_PazyFFWT(hingeNode=hingeNode,hingeAngle=hingeAngle,flareAngle=flareAngle,kSpring=kSpring,airspeed=U,p0=[0;0;θ],g=g)
# plt = plot_undeformed_assembly(pazyFFWT)
# display(plt)

# System solver
σ0 = 1.0
maxIter = 50
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter)

# Create and solve problem
problem = create_SteadyProblem(model=pazyFFWT,systemSolver=NR)
solve!(problem)

# Get outputs
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1_mainWing = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:hingeNode-1]...)
x1_wingTip = x1_mainWing[end] .+ vcat([vcat(problem.model.beams[2].elements[e].x1_n1,problem.model.beams[2].elements[e].x1_n2) for e in 1:16-hingeNode]...)
x1 = vcat(x1_mainWing,x1_wingTip)
u1_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...)
u3_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
p2_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
M2_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)

# Plots
# ------------------------------------------------------------------------------
# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath="/test/outputs/figures/PazyFFWTsteady/PazyFFWTsteady_deformation.pdf")
display(deformationPlot)
# OOP displacement
plt1 = plot((r_n1.+u1_of_x1)/x1[end], (r_n3.+u3_of_x1)/x1[end], aspect_ratio=:equal, lw=2, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyFFWTsteady/PazyFFWTsteady_u3.pdf"))
# Bending angle
plt2 = plot(x1/x1[end], p2_of_x1, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyFFWTsteady/PazyFFWTsteady_p2.pdf"))
# Bending moment
plt3 = plot(x1/x1[end], M2_of_x1, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/PazyFFWTsteady/PazyFFWTsteady_M2.pdf"))

println("Finished PazyFFWTsteady.jl")