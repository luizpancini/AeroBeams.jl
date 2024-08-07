using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = 12.0
EA = 30e6
EIy = 30e6*1/12
∞ = 1e14
stiffnessMatrix = diagm([EA,∞,∞,∞,EIy,∞])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipMoment = create_BC(name="tipMoment",beam=beam,node=nElem+1,types=["M2A"],values=[2*π*EIy/L])

# Model
tipMomentCantilever = create_Model(name="tipMomentCantilever",beams=[beam],BCs=[clamp,tipMoment],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tipMomentCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][nElem].θ_n2 for i in 1:length(σVector)]

# Plot internal bending moment
plot_steady_outputs(problem,outputs=["M2"],save=true,saveFolder="/test/outputs/figures/tipMomentCantilever/")

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,save=true,savePath="/test/outputs/figures/tipMomentCantilever/tipMomentCantilever_deformation.pdf")
display(deformationPlot)

# Plot normalized displacements over load steps
gr()
y = [1.0 .+ tip_u1/L, tip_u3/L, tip_angle/π]
labels = ["\$1+u_1/L\$" "\$u_3/L\$" "\$\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
plt1 = plot(xlabel="\$ML/(2\\pi EI)\$", ylabel="\$1+u_1/L, u_3/L, \\theta/L\$", title="Tip generalized displacements")
plot!(σVector, y, palette=colors, lw=2, label=false)
halfNσ = round(Int,length(σVector)/2)
for i=1:3
    if i==1
        annotate!(σVector[halfNσ], y[i][halfNσ], text(labels[i], :bottom, :left, colors[i]))
    else
        annotate!(σVector[halfNσ], y[i][halfNσ], text(labels[i], :bottom, :right, colors[i]))
    end
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/tipMomentCantilever/tipMomentCantilever_summary.pdf"))

println("Finished tipMomentCantilever.jl")