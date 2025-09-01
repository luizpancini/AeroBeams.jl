using Plots, ColorSchemes

# Run the script
include("../examples/S10PazyFlutterOnsetPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/S10PazyFlutterOnsetPitchRange"
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
plt = plot(xlabel="Tip displacement [m]", ylabel="Flutter speed [m/s]", xlims=[0,0.3+0.01], ylims=[0,80], xticks=vcat(0:0.05:0.3), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:], flutterOnsetSpeed[c,:], c=configColors[c], lw=lw, label=false)
end
display(plt)
savefig(plt,string(absPath,"/S10PazyFlutterOnsetPitchRange_Uf_vs_tipDisp.pdf"))

println("Finished S10PazyFlutterOnsetPitchRangePlotGenerator.jl")