using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

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

# Load reference solution
u1Ref = readdlm(string(pwd(),"/test/referenceData/distributedLoadCantilever/u1.txt"))
u3Ref = readdlm(string(pwd(),"/test/referenceData/distributedLoadCantilever/u3.txt"))
θRef = readdlm(string(pwd(),"/test/referenceData/distributedLoadCantilever/theta.txt"))

# Plot deformed state
relPath = "/test/outputs/figures/distributedLoadCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/distributedLoadCantilever_deformation.pdf"))
display(deformationPlot)

# Plot normalized displacements over load steps
lw = 2
ms = 4
x = [-tip_u1/L, tip_u3/L, -tip_angle/π]
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()
plt1 = plot(xlabel="\$-u_1/L, u_3/L, -\\theta/L\$", ylabel="\$q [kN]\$", title="Tip generalized displacements",legend=:bottomright)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*q, palette=colors, lw=lw, label=false)
scatter!([u1Ref[1,:],u3Ref[1,:],θRef[1,:]], [u1Ref[2,:],u3Ref[2,:],θRef[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/distributedLoadCantilever_disp.pdf"))

println("Finished distributedLoadCantilever.jl")