using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyEigenStructuralSweepCorrectionDispRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyEigenStructuralSweepCorrectionDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad([:blue,:red,:gold], nModes, categorical=true)
ls = [:solid,:dash,:dashdot]
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
gr()

# Frequency vs. tip displacement
plt = plot(xlabel="Tip displacement [% semispan]", ylabel="Frequency [Hz]", xlims=[0,50], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,Λ) in enumerate(ΛRange)
    Λstr = string(round(Int,Λ*180/pi))
    plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[i], label="\$\\Lambda=$Λstr ^\\circ\$")
end
for (i,Λ) in enumerate(ΛRange)
    for mode in 1:nModes
        plot!(tipOOP[i,:]/L*100, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, ls=ls[i], label=false)
    end
end
display(plt)
savefig(string(absPath,"/sweptPazyEigenStructuralSweepCorrectionDispRange.pdf"))

println("Finished sweptPazyEigenStructuralSweepCorrectionDispRange.jl")