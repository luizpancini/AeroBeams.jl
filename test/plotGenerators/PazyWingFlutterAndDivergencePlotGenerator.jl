using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingFlutterAndDivergence.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingFlutterAndDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

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
savefig(string(absPath,"/PazyWingFlutterAndDivergence_Vgf.pdf"))

println("Finished PazyWingFlutterAndDivergencePlotGenerator.jl")