using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutterPrecurvatureRange.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutterPrecurvatureRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(kRange)))
lw = 2
ms = 3
gr()

# Flutter speed vs. tip load
plt1 = plot(xlabel="Tip load [N]", ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,40])
for (i,k) in enumerate(kRange)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw=lw, label="AeroBeams - \$k\$=$k")
end
plot!(flutterSpeedVsTipLoadk0[1,:], flutterSpeedVsTipLoadk0[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plot!(flutterSpeedVsTipLoadk2[1,:], flutterSpeedVsTipLoadk2[2,:], c=:black, ls=:dash, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_flutter.pdf"))

println("Finished SMWFlutterPrecurvatureRangePlotGenerator.jl")