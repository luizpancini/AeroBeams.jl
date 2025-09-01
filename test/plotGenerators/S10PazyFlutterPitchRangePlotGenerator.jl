using Plots, ColorSchemes

# Run the script
include("../examples/S10PazyFlutterPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/S10PazyFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# V-g-f
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,105], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
for (c,config) in enumerate(tipMassConfigs)
    for (i,θ) in enumerate(θRange)
        for mode in 1:nModes
            plot!(URange, modeFrequencies[c,i,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        end
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for (i,θ) in enumerate(θRange)
        for mode in 1:nModes
            plot!(URange, modeDampingRatios[c,i,mode], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        end
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vf)
    display(plt_Vg)
    display(plt_Vgf)
    savefig(plt_Vf,string(absPath,"/S10PazyFlutterPitchRange_freq_",config,".pdf"))
    savefig(plt_Vg,string(absPath,"/S10PazyFlutterPitchRange_damp_",config,".pdf"))
    savefig(plt_Vgf,string(absPath,"/S10PazyFlutterPitchRange_Vgf_",config,".pdf"))
end

println("Finished S10PazyFlutterPitchRange.jl")