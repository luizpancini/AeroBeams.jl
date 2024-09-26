using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingFlutter.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes
modesPlot = plot_mode_shapes(problem[end],scale=0.5,legendPos=(0.25,0.2),view=(30,30),save=true,savePath=string(relPath,"/PazyWingFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

# Normalized deformed wingspan
plt0 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt0)
savefig(string(absPath,"/PazyWingFlutter_disp.pdf"))

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingFlutter_Vgf.pdf"))

# Frequencies and dampings vs tip OOP displacement
plt21 = plot(ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(tip_OOP/L*100, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw,  label=false)
end
plt22 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Damping Ratio", ylims=[-0.15,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(tip_OOP/L*100, modeDampingRatios[mode], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt2 = plot(plt21,plt22, layout=(2,1))
display(plt2)
savefig(string(absPath,"/PazyWingFlutter_OOPgf.pdf"))

println("Finished PazyWingFlutterPlotGenerator.jl")