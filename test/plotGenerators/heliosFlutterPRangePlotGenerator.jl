using Plots, ColorSchemes

# Run the script
include("../examples/heliosFlutterPRange.jl")

# Set paths
relPath = "/test/outputs/figures/heliosFlutterPRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes
modesPlot = plot_mode_shapes(eigenProblem,scale=10,view=(30,30),legendPos=:outerright,nModes=6,save=true,savePath=string(relPath,"/heliosFlutterPRange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = get(colorschemes[:jet1], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(PRange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt12 = plot(xlabel="Payload [lb]", ylabel="Damping [1/s]", ylims=[-5.0,1.0], legend=:bottomleft)
for mode in 1:nModes
    scatter!(PRange, modeDampings[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosFlutterPRange_Vgf.pdf"))

# Root locus 
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosFlutterPRange_rootlocus.pdf"))

# Root locus (zoom)
plt3 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt3)
savefig(string(absPath,"/heliosFlutterPRange_rootlocuszoom.pdf"))

# Root locus (phugoid zoom)
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.1,0.2], ylims=[0,0.6])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/heliosFlutterPRange_Phzoom.pdf"))

println("Finished heliosFlutterPRangePlotGenerator.jl")