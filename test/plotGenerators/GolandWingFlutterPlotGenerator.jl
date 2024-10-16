using Plots, ColorSchemes

# Run the script
include("../examples/GolandWingFlutter.jl")

# Set paths
relPath = "/test/outputs/figures/GolandWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 2
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.25,0.1],legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/GolandWingFlutter_Vgf.pdf"))

println("Finished GolandWingFlutterPlotGenerator.jl")