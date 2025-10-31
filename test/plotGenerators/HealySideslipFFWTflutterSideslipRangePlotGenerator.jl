using Plots, ColorSchemes

# Run the script
include("../examples/HealySideslipFFWTflutterSideslipRange.jl")

# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTflutterSideslipRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
θColors = cgrad(:rainbow, length(θRange), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 12
fs = 14
lfs = 12
lw = 2
ms = 3
gr()
θLabels = ["\$\\theta=$(round(Int,θ*180/π))^\\circ \$" for θ in θRange]

# V-g-f for all configurations
plt_Vf = plot(ylabel="Frequency [Hz]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, xlims=180/π.*extrema(βRange), ylims=[0,20], legend=:topright)
plt_Vg = plot(xlabel="Sideslip angle [deg]", ylabel="Damping Ratio", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, xlims=180/π.*extrema(βRange), ylims=[-0.25,0.2], legend=:topleft)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        currentLabel = mode == 1 ? θLabels[i] : false
        plot!(plt_Vf,βRange*180/π, modeFrequencies[i,mode]/(2*π), c=θColors[i], lw=lw, label=currentLabel)
    end
    for mode in 1:nModes
        plot!(plt_Vg,βRange*180/π, modeDampingRatios[i,mode], c=θColors[i], lw=lw, label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/HealySideslipFFWTflutterSideslipRange.pdf"))

println("Finished HealySideslipFFWTflutterSideslipRangePlotGenerator.jl")