# Run the script
include("../examples/tipMomentCantilever.jl")

# Plot deformed state
relPath = "/test/outputs/figures/tipMomentCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/tipMomentCantilever_deformation.pdf"))
display(deformationPlot)

# Plot internal bending moment and strain
gr()
plot_steady_outputs(problem,outputs=["M2","κ2"],save=true,saveFolder=string(relPath,"/"))

# Plot normalized displacements over load steps
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
savefig(string(absPath,"/tipMomentCantilever_disp.pdf"))

println("Finished tipMomentCantileverPlotGenerator.jl")