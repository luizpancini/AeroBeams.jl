using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTlockedFlutterURange.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTlockedFlutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at lowest airspeed
modesPlot = plot_mode_shapes(problem[1],scale=1,frequencyLabel="frequency",view=(30,30),modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/HealyBaselineFFWTlockedFlutterURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]", xlims=[0,40], ylims=[0,80])
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,40], ylims=[-0.1,0.05], legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/HealyBaselineFFWTlockedFlutterURange.pdf"))

println("Finished HealyBaselineFFWTlockedFlutterURangePlotGenerator.jl")