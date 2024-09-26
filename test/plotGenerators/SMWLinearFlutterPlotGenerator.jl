using Plots, ColorSchemes

# Run the script
include("../examples/SMWLinearFlutter.jl")

# Set paths
relPath = "/test/outputs/figures/SMWLinearFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem[end],scale=5,view=(30,30),legendPos=:best,frequencyLabel="frequency",save=true,savePath=string(relPath,"/SMWLinearFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 4
msw = 0
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], mc=modeColors[mode], ms=ms, msw=msw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.2,0.1])
for mode in 1:nModes
    scatter!(URange, modeDampingRatios[mode], mc=modeColors[mode], ms=ms, msw=msw, label=false)
    scatter!([NaN], [NaN], mc=modeColors[mode], ms=ms, msw=msw, label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWLinearFlutter_Vgf.pdf"))

println("Finished SMWLinearFlutterPlotGenerator.jl")