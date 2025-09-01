using Plots, ColorSchemes

# Run the script
include("../examples/S20PazyFlutterOnsetPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/S20PazyFlutterOnsetPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
configColors = cgrad(:rainbow, length(tipMassConfigs), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 1
gr()

# Flutter onset speed vs. tip displacement
plt = plot(xlabel="Tip displacement [m]", ylabel="Flutter speed [m/s]", xlims=[0,0.25+0.01], ylims=[0,80], xticks=vcat(0:0.05:0.25), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:black, label="Exp. (gyros)")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:square, ms=ms, msw=msw, mc=:white, msc=:black, label="Exp. (cameras)")
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:], flutterOnsetSpeed[c,:], c=configColors[c], lw=lw, label=false)
    if config == "LE"
        plot!(Uf_vs_tipdisp_LE_g[1,:], Uf_vs_tipdisp_LE_g[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
        plot!(Uf_vs_tipdisp_LE_c[1,:], Uf_vs_tipdisp_LE_c[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:square, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == "TE"
        plot!(Uf_vs_tipdisp_TE_g[1,:], Uf_vs_tipdisp_TE_g[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
        plot!(Uf_vs_tipdisp_TE_c[1,:], Uf_vs_tipdisp_TE_c[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:square, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end
display(plt)
savefig(plt,string(absPath,"/S20PazyFlutterOnsetPitchRange_Uf_vs_tipDisp.pdf"))

println("Finished S20PazyFlutterOnsetPitchRangePlotGenerator.jl")