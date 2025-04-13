using Plots, ColorSchemes

# Run the script
include("../examples/typicalSectionFlutterAndDivergence.jl")

# Set paths
relPath = "/test/outputs/figures/typicalSectionFlutterAndDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes+1))
modeLabels = ["Plunge" "Pitch"]
lw = 3
ms = 2
msw = 0
gr()

# V-g-f
plt11 = plot(ylabel="Frequency ratio", xlims=[0,URange[end]/(b*ωα)], ylims=[0,1.02])
for mode in 1:nModes
    plot!(URange/(b*ωα), modeFrequencies[mode]/ωα, c=modeColors[mode], lw=lw,  label=false)
end
plot!(freqsRef[1,:], freqsRef[2,:], ls=:dash, lw=1, marker=:circle, ms=2, msw=msw, c=modeColors[end], label=false)
plt12 = plot(xlabel="Normalized airspeed, \$U/b\\omega_\\alpha\$", ylabel="Damping Ratio", xlims=[0,URange[end]/(b*ωα)], ylims=[-0.15,0.15],legend=:topleft)
for i in eachindex(URange)
    for j in eachindex(dampingsNonOscillatory[i])
        scatter!([URange[i]/(b*ωα)], [dampingsNonOscillatory[i][j]/ωα], c=:black, ms=ms, msw=msw, label=false)
    end
end
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Non-oscillatory")
for mode in 1:nModes
    plot!(URange/(b*ωα), modeDampings[mode]/ωα, c=modeColors[mode], lw=lw, label=modeLabels[mode])
end
plot!(dampsRef[1,:], dampsRef[2,:], ls=:dash, lw=1, marker=:circle, ms=2, msw=msw, c=modeColors[end], label="Hodges & Pierce (2011) - \$p\$ method")
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/typicalSectionFlutter_Vgf.pdf"))

println("Finished typicalSectionFlutterAndDivergencePlotGenerator.jl")