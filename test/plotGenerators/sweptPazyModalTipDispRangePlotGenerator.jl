using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyModalTipDispRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyModalTipDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lfs = 12
lw = 2
ls = [:solid :dash :dot]
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Normalized tip displacement", ylabel="Frequency [Hz]", xlims=[0,0.5], ylims=[0,50], legendfontsize=lfs, colorbar_title="Mode")
for (i,Λ) in enumerate(ΛRange)
    plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[i], label="\$\\Lambda=\$$(round.(Int,Λ*180/π))\$^\\circ\$")
end
for (i,Λ) in enumerate(ΛRange)
    for mode in 1:nModes
        plot!(tipOOP[i,:]/L, modeFrequencies[i,mode], c=colors, lw=lw, ls=ls[i], lz=mode, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/sweptPazyModalTipDispRange.pdf"))

println("Finished sweptPazyModalTipDispRange.jl")