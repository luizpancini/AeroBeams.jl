using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTfreeFlutterFlareRangeURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
flareColors = cgrad(:rainbow, length(ΛRange), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 4
msw = 0
gr()

# V-g-fs
for (i,Λ) in enumerate(ΛRange)
    Λdeg = round(Int,Λ*180/π)
    plt11 = plot(ylabel="Frequency [Hz]", xlims=[0,40], ylims=[0,25],title=string("\$\\Lambda = ",Λdeg,"\\degree\$"))
    for mode in 1:nModes
        scatter!(URange, modeFrequencies[i,mode]/(2*π), c=modeColors[mode], ms=ms, msw=msw, label=false)
    end
    plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,40], ylims=[-0.5,0.5], legend=:topleft)
    for mode in 1:nModes
        scatter!(URange, modeDampingRatios[i,mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
    end
    plt1 = plot(plt11,plt12, layout=(2,1))
    display(plt1)
    savefig(string(absPath,"/HealyBaselineFFWTfreeFlutterFlareRangeURange_Vgf_Lambda_",Λdeg,".pdf"))
end

# V-g-f - zoom on first two modes
plt21 = plot(ylabel="Frequency [Hz]", xlims=[0,40], ylims=[0,5], legend=:bottomright)
plt22 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,40], ylims=[-1,0.5])
plot!(plt21, [NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
scatter!(plt21, [NaN],[NaN], c=:black, ms=ms, msw=msw, label="AeroBeams")
for (i,Λ) in enumerate(ΛRange)
    if i==1
        mode1_damp = flare10_mode1_damp
        mode2_damp = flare10_mode2_damp
        mode1_freq = flare10_mode1_freq
        mode2_freq = flare10_mode2_freq
    elseif i==2
        mode1_damp = flare15_mode1_damp
        mode2_damp = flare15_mode2_damp
        mode1_freq = flare15_mode1_freq
        mode2_freq = flare15_mode2_freq
    else   
        mode1_damp = flare20_mode1_damp
        mode2_damp = flare20_mode2_damp
        mode1_freq = flare20_mode1_freq
        mode2_freq = flare20_mode2_freq 
    end
    for mode in [1,2]
        scatter!(plt21, URange, modeFrequencies[i,mode]/(2*π), c=flareColors[i], ms=ms, msw=msw, label=false)
        scatter!(plt22, URange, modeDampingRatios[i,mode], c=flareColors[i], ms=ms, msw=msw, label=false)
    end
    plot!(plt21, mode1_freq[1,:], mode1_freq[2,:], c=flareColors[i], ls=:dash, lw=lw, label=false)
    plot!(plt21, mode2_freq[1,:], mode2_freq[2,:], c=flareColors[i], ls=:dash, lw=lw, label=false)
    plot!(plt22, mode1_damp[1,:], mode1_damp[2,:], c=flareColors[i], ls=:dash, lw=lw, label=false)
    plot!(plt22, mode2_damp[1,:], mode2_damp[2,:], c=flareColors[i], ls=:dash, lw=lw, label=false)
    plot!(plt21,[NaN],[NaN], c=flareColors[i], lw=lw, label="\$\\Lambda = $(round(Int,Λ*180/π)) \\degree\$")
end
plt2 = plot(plt21,plt22, layout=(2,1))
display(plt2)
savefig(string(absPath,"/HealyBaselineFFWTfreeFlutterFlareRangeURange_Vgfzoom.pdf"))

println("Finished HealyBaselineFFWTfreeFlutterFlareRangeURangePlotGenerator.jl")