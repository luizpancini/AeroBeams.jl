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
lw = 2
ms = 2
gr()

# V-g-f
plt11 = plot(ylabel="Frequency ratio", ylims=[0,1.01])
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/ωα, c=modeColors[mode], lw=lw,  label=false)
end
for r in 2:axes(freqsRef, 1)[end]
    plot!(freqsRef[1,:], freqsRef[r,:], ls=:dash, c=modeColors[end], lw=1, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.15],legend=:topleft)
for i in eachindex(URange)
    for j in eachindex(dampingsNonOscillatory[i])
        scatter!([URange[i]], [dampingsNonOscillatory[i][j]/ωα], c=:black, ms=ms, msw=0, label=false)
    end
end
scatter!([NaN], [NaN], c=:black, ms=ms, msw=0, label="Non-oscillatory")
for mode in 1:nModes
    plot!(URange, modeDampings[mode]/ωα, c=modeColors[mode], lw=lw, label=modeLabels[mode])
end
for r in 2:axes(dampsRef, 1)[end]
    plot!(dampsRef[1,:], dampsRef[r,:], ls=:dash, c=modeColors[end], lw=1, label=false)
end
plot!([NaN], [NaN], ls=:dash, c=modeColors[end], lw=1, label="\$p\$ method")
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/typicalSectionFlutter_Vgf.pdf"))

println("Finished typicalSectionFlutterAndDivergencePlotGenerator.jl")