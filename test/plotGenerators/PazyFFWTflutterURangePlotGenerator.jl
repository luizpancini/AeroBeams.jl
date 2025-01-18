using Plots, ColorSchemes

# Run the script
include("../examples/PazyFFWTflutterURange.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTflutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at lowest airspeed
modesPlot = plot_mode_shapes(problem[1],scale=2,view=(30,30),modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/PazyFFWTflutterURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]", xlims=[0,URange[end]+1e-3], ylims=[0,50])
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,URange[end]+1e-3], ylims=[-0.1,0.05], legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyFFWTflutterURange_Vgf.pdf"))

# Coast (or fold) angle
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Hinge angle [deg]", xlims=[0,URange[end]+1e-3], ylims=[-135,135], yticks=-135:45:135)
plot!(URange, -ϕHinge, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyFFWTflutterURange_phi.pdf"))

println("Finished PazyFFWTflutterURangePlotGenerator.jl")