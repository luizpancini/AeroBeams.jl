using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyFlutterPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
tipLossStr = ifelse(hasTipCorrection, tipLossType, "none")
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# Assumed experimental air density
ρref = 1.225

# V-g-f
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,100], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
scatter!([NaN], [NaN], ms=ms, msw=msw, mc=:white, msc=:black, label="Exp.")
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        if θ == 0 && mode <= size(freqs_aoa0_ref,1)-1
            scatter!(sqrt.(2*freqs_aoa0_ref[1,:]/ρref), freqs_aoa0_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 1*π/180 && mode <= size(freqs_aoa1_ref,1)-1
            scatter!(sqrt.(2*freqs_aoa1_ref[1,:]/ρref), freqs_aoa1_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 3*π/180 && mode <= size(freqs_aoa3_ref,1)-1
            scatter!(sqrt.(2*freqs_aoa3_ref[1,:]/ρref), freqs_aoa3_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 5*π/180 && mode <= size(freqs_aoa5_ref,1)-1
            scatter!(sqrt.(2*freqs_aoa5_ref[1,:]/ρref), freqs_aoa5_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 7*π/180 && mode <= size(freqs_aoa7_ref,1)-1
            scatter!(sqrt.(2*freqs_aoa7_ref[1,:]/ρref), freqs_aoa7_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 10*π/180 && mode <= size(freqs_aoa10_ref,1)-1
            scatter!(sqrt.(2*freqs_aoa10_ref[1,:]/ρref), freqs_aoa10_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)      
        end
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[i,mode], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vf)
display(plt_Vg)
display(plt_Vgf)
savefig(plt_Vf,string(absPath,"/sweptPazyFlutterPitchRange_freq_",tipLossStr,"_Lambda",string(round(Int,Λ*180/pi)),".pdf"))
savefig(plt_Vg,string(absPath,"/sweptPazyFlutterPitchRange_damp_",tipLossStr,"_Lambda",string(round(Int,Λ*180/pi)),".pdf"))
savefig(plt_Vgf,string(absPath,"/sweptPazyFlutterPitchRange_Vgf_",tipLossStr,"_Lambda",string(round(Int,Λ*180/pi)),".pdf"))

println("Finished sweptPazyFlutterPitchRange.jl")