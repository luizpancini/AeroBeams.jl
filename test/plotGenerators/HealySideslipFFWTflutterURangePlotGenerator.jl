using Plots, ColorSchemes

# Run the script
include("../examples/HealySideslipFFWTflutterURange.jl")

# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTflutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at lowest airspeed
modesPlot = plot_mode_shapes(problem[1],scale=2,view=(30,30),modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/HealySideslipFFWTflutterURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 12
fs = 14
lfs = 12
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]", xlims=[URange[1],URange[end]+1e-3], ylims=[0,40], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]+1e-3], ylims=[-0.1,0.05], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/HealySideslipFFWTflutterURange.pdf"))

println("Finished HealySideslipFFWTflutterURangePlotGenerator.jl")