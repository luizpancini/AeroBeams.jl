using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutterTipDispRange.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutterTipDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 3
gr()

# Flutter speed and frequency vs. tip displacement
plt11 = plot(ylabel="Flutter speed [m/s]", xlims=[0,3], ylims=[0,35], legend=:bottomleft)
plot!(flutterTipDisp, flutterSpeed, c=:black, lw=lw, label="AeroBeams")
plot!(flutterSpeedRef[1,:], flutterSpeedRef[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plt12 = plot(xlabel="Tip displacement [m]", ylabel="Flutter frequency [rad/s]", xlims=[0,3], ylims=[0,35])
plot!(flutterTipDisp, flutterFreq, c=:black, lw=lw, label=false)
plot!(flutterFreqRef[1,:], flutterFreqRef[2,:], c=:black, ls=:dash, lw=lw, label=false)
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWFlutterTipDispRange_flutterVsDisp.pdf"))

# Complete flutter speed curve
plt2 = plot(xlabel="Tip displacement [m]", ylabel="Flutter speed [m/s]", xlims=[-3,3], ylims=[0,35], legend=:bottomleft)
plot!(flutterTipDisp, flutterSpeed, c=:black, lw=lw, label="AeroBeams")
plot!(flutterSpeedVsDispRef[1,:], flutterSpeedVsDispRef[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
display(plt2)
savefig(string(absPath,"/SMWFlutterTipDispRange_flutterCurve.pdf"))

println("Finished SMWFlutterTipDispRangePlotGenerator.jl")