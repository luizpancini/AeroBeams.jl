using Plots, ColorSchemes

# Run the script
include("../examples/heliosWingFlutter.jl")

# Set paths
relPath = "/test/outputs/figures/heliosWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes
modesPlot = plot_mode_shapes(eigenProblem,scale=5,view=(30,30),frequencyLabel="frequency",save=true,savePath=string(relPath,"/heliosWingFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = get(colorschemes[:jet1], LinRange(0, 1, nModes))
lw = 2
ms = 3
msw = 0
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]", ylims=[0,15])
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-2.0,0.25], legend=:bottomleft)
for mode in 1:nModes
    scatter!(URange, modeDampingRatios[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosWingFlutter_Vgf.pdf"))

# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-2.0,0.25] ,ylims=[0,15])
for mode in 1:nModes
    scatter!(modeDampingRatios[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosWingFlutter_rootlocus.pdf"))

println("Finished heliosWingFlutterPlotGenerator.jl")