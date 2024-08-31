using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingFlutterTipMassRange.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingFlutterTipMassRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
lstyles = [:solid :dash :dot]
gr()

# V-g-f for all configurations
plt11 = plot(ylabel="Frequency [Hz]")
for c in configurations
    for mode in 1:nModes
        plot!(URange, modeFrequencies[c,mode]/(2*π), c=modeColors[mode], ls=lstyles[c], lw=lw, label=false)
    end
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.1,0.05], legend=:topleft)
plot!([NaN], [NaN], c=:black, ls=lstyles[1], lw=lw, label="Clean wing")
plot!([NaN], [NaN], c=:black, ls=lstyles[2], lw=lw, label="LE weight")
plot!([NaN], [NaN], c=:black, ls=lstyles[3], lw=lw, label="TE weight")
for c in configurations
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[c,mode], c=modeColors[mode], ls=lstyles[c], lw=lw, label=false)
    end
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingFlutterTipMassRange.pdf"))

println("Finished PazyWingFlutterTipMassRangePlotGenerator.jl")