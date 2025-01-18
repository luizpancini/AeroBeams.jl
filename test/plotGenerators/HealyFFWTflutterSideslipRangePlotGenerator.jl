using Plots, ColorSchemes

# Run the script
include("../examples/HealyFFWTflutterSideslipRange.jl")

# Set paths
relPath = "/test/outputs/figures/HealyFFWTflutterSideslipRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f for all configurations
plt11 = plot(ylabel="Frequency [Hz]", xlims=180/π*[βRange[1],βRange[end]+1e-3], ylims=[0,20])
for mode in 1:nModes
    plot!(βRange*180/π, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Sideslip angle [deg]", ylabel="Damping Ratio", xlims=180/π*[βRange[1],βRange[end]+1e-3], ylims=[-0.2,0.1], legend=:topleft)
for mode in 1:nModes
    plot!(βRange*180/π, modeDampingRatios[mode], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/HealyFFWTflutterSideslipRange.pdf"))

println("Finished HealyFFWTflutterSideslipRangePlotGenerator.jl")