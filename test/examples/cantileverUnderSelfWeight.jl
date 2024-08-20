using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L,b,H = 1,1e-2,1e-3
E,ρ = 200e9,7.8e3
A,Iy = b*H,b*H^3/12
∞ = 1e12
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*Iy,∞])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,0,0,0])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
cantileverUnderSelfWeight = create_Model(name="cantileverUnderSelfWeight",beams=[beam],BCs=[clamp],gravityVector=[0,0,-9.81])

# Set system solver options
σ0 = 0
σstep = 0.05
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=cantileverUnderSelfWeight,systemSolver=NR)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Plots
# ------------------------------------------------------------------------------
relPath = "/test/outputs/figures/cantileverUnderSelfWeight"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/cantileverUnderSelfWeight_deformation.pdf"))
display(deformationPlot)
# u1
gr()
plt0 = plot(x1/L, u1/L, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$u_1/L\$")
display(plt0)
savefig(string(absPath,"/cantileverUnderSelfWeight_u1.pdf"))
# u3
plt1 = plot(x1/L, u3/L, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/cantileverUnderSelfWeight_u3.pdf"))
# F3
plt2 = plot(x1/L, F3, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/cantileverUnderSelfWeight_F3.pdf"))
# M2
plt3 = plot(x1/L, M2, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/cantileverUnderSelfWeight_M2.pdf"))

println("Finished cantileverUnderSelfWeight.jl")