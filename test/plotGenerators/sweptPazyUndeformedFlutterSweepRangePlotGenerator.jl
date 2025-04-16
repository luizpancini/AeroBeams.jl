using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyUndeformedFlutterSweepRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyUndeformedFlutterSweepRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad([:blue,:red,:green], nModes, categorical=true)
ts = 10
fs = 15
lfs = 8
lw = 2
ms = 4
gr()

# V-g-f
for (i,Λ) in enumerate(ΛRange)
    Λstr = string(round(Int,Λ*180/pi))
    ULim = i == 1 ? 121 : 90
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[10,ULim], ylims=[0,50], xticks=0:10:ULim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legendposition=(0.75,0.95))
    if i==1
        plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams (VLM)")
        scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Nastran (Technion)")
    end
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, label=false)
    end
    if i == 1
        scatter!(undef_freq_Lambda0_NastranTechnion[1,1:46],undef_freq_Lambda0_NastranTechnion[2,1:46], c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda0_NastranTechnion[1,48:93],undef_freq_Lambda0_NastranTechnion[2,48:93], c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda0_NastranTechnion[1,95:end],undef_freq_Lambda0_NastranTechnion[2,95:end], c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 3
        scatter!(undef_freq_Lambda20_NastranTechnion[1,1:31],undef_freq_Lambda20_NastranTechnion[2,1:31], c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda20_NastranTechnion[1,33:63],undef_freq_Lambda20_NastranTechnion[2,33:63], c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda20_NastranTechnion[1,65:end],undef_freq_Lambda20_NastranTechnion[2,65:end], c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 4
        scatter!(undef_freq_Lambda30_NastranTechnion[1,1:31],undef_freq_Lambda30_NastranTechnion[2,1:31], c=modeColors[1], ms=ms, msw=msw, label=false)    
        scatter!(undef_freq_Lambda30_NastranTechnion[1,33:63],undef_freq_Lambda30_NastranTechnion[2,33:63], c=modeColors[2], ms=ms, msw=msw, label=false)  
        scatter!(undef_freq_Lambda30_NastranTechnion[1,65:end],undef_freq_Lambda30_NastranTechnion[2,65:end], c=modeColors[3], ms=ms, msw=msw, label=false)  
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[10,ULim], ylims=[-0.1,0.075], xticks=0:10:ULim, yticks=-0.1:0.05:0.1, tickfont=font(ts), guidefont=font(fs))
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[i,mode], c=modeColors[mode], lw=lw, label=false)
    end
    if i == 1
        scatter!(undef_damp_Lambda0_NastranTechnion[1,1:6],undef_damp_Lambda0_NastranTechnion[2,1:6]/2, c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda0_NastranTechnion[1,8:37],undef_damp_Lambda0_NastranTechnion[2,8:37]/2, c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda0_NastranTechnion[1,39:end],undef_damp_Lambda0_NastranTechnion[2,39:end]/2, c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 3
        scatter!(undef_damp_Lambda20_NastranTechnion[1,1:4],undef_damp_Lambda20_NastranTechnion[2,1:4]/2, c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda20_NastranTechnion[1,6:33],undef_damp_Lambda20_NastranTechnion[2,6:33]/2, c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda20_NastranTechnion[1,35:end],undef_damp_Lambda20_NastranTechnion[2,35:end]/2, c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 4
        scatter!(undef_damp_Lambda30_NastranTechnion[1,1:5],undef_damp_Lambda30_NastranTechnion[2,1:5]/2, c=modeColors[1], ms=ms, msw=msw, label=false)    
        scatter!(undef_damp_Lambda30_NastranTechnion[1,6:36],undef_damp_Lambda30_NastranTechnion[2,6:36]/2, c=modeColors[2], ms=ms, msw=msw, label=false) 
        scatter!(undef_damp_Lambda30_NastranTechnion[1,38:end],undef_damp_Lambda30_NastranTechnion[2,38:end]/2, c=modeColors[3], ms=ms, msw=msw, label=false) 
    end
    annotate!([30],[0.03], text("\$\\Lambda=$Λstr ^\\circ\$", 30, :black))
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(string(absPath,"/sweptPazyUndeformedFlutterSweepRange_Vgf_Lambda",Λstr,".pdf"))
end

println("Finished sweptPazyUndeformedFlutterSweepRangePlotGenerator.jl")