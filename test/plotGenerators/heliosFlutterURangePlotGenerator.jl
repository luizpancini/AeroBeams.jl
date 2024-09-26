using Plots, ColorSchemes

# Run the script
include("../examples/heliosFlutterURange.jl")

# Set paths
relPath = "/test/outputs/figures/heliosFlutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes
modesPlot = plot_mode_shapes(eigenProblem,scale=5,view=(30,30),save=true,savePath=string(relPath,"/heliosFlutterURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = get(colorschemes[:jet1], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-1.0,0.25], legend=:bottomleft)
for mode in 1:nModes
    scatter!(URange, modeDampingRatios[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosFlutterURange_Vgf.pdf"))

# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-1.0,0.25])
for mode in 1:nModes
    scatter!(modeDampingRatios[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosFlutterURange_rootlocus.pdf"))

println("Finished heliosFlutterURangePlotGenerator.jl")