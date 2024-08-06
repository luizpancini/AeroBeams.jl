using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = 1
E = 210e6
A,Iy = 20e-4,5/3*1e-8
∞ = 1e12
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*Iy,∞])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
q = 100
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
add_loads_to_beam!(beam,loadTypes=["ff_b_of_x1t"],loadFuns=[(x1,t)->[0; 0; q]])

# Model
distributedLoadCantilever = create_Model(name="distributedLoadCantilever",beams=[beam],BCs=[clamp])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=distributedLoadCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][nElem].θ_n2 for i in 1:length(σVector)]

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,save=true,savePath="/test/outputs/figures/distributedLoadCantilever/distributedLoadCantilever_deformation.pdf")
display(deformationPlot)

# Plot normalized displacements over load steps
x = [-tip_u1/L, tip_u3/L, -tip_angle/π]
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()
plt1 = plot(xlabel="\$-u_1/L, u_3/L, -\\theta/L\$", ylabel="\$q [kN]\$", title="Tip generalized displacements")
plot!(x, σVector*q, palette=colors, lw=2, label=false)
halfNσ = round(Int,length(σVector)/2)
tqNσ = round(Int,length(σVector)*3/4)
annotate!(x[1][halfNσ], σVector[halfNσ]*q, text(labels[1], :top, :left, colors[1]))
annotate!(x[2][tqNσ], σVector[tqNσ]*q, text(labels[2], :top, :right, colors[2]))
annotate!(x[3][tqNσ], σVector[tqNσ]*q, text(labels[3], :top, :left, colors[3]))
display(plt1)

println("Finished distributedLoadCantilever.jl")