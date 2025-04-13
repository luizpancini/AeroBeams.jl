using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyModalTipDispRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyModalTipDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ls = [:solid :dash :dot]
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Tip OOP disp. [% semispan]", ylabel="Frequency [Hz]", xlims=[0,50], ylims=[0,50])
for (i,Λ) in enumerate(ΛRange)
    for mode in 1:nModes
        plot!(tipOOP[i,:]/L*100, modeFrequencies[i,mode], c=colors[mode], lw=lw, ls=ls[i], label=false)
    end
end
display(plt1)
savefig(string(absPath,"/sweptPazyModalTipDispRange.pdf"))

println("Finished sweptPazyModalTipDispRange.jl")