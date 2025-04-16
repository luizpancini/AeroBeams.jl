using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyFlutterPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
tipLossStr = ifelse(hasTipCorrection, tipLossType, "none")
modeColors = cgrad([:blue,:red,:green], nModes, categorical=true)
ts = 10
fs = 16
lw = 2
ms = 3
gr()

# V-g-f
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,50], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], tickfont=font(ts), guidefont=font(fs), legend=:topleft)
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[i,mode], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vg)
display(plt_Vgf)
savefig(plt_Vg,string(absPath,"/sweptPazyFlutterPitchRange_damp_",tipLossStr,"_Lambda",string(round(Int,Λ*180/pi)),".pdf"))
savefig(plt_Vgf,string(absPath,"/sweptPazyFlutterPitchRange_Vgf_",tipLossStr,"_Lambda",string(round(Int,Λ*180/pi)),".pdf"))

println("Finished sweptPazyFlutterPitchRange.jl")