using Plots, ColorSchemes

# Run the script
include("../examples/sweptWingBendingDivergenceSweepRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptWingBendingDivergenceSweepRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
gr()

# Divergence speed vs. sweep angle
plt_UD = plot(xlabel="Sweep angle [deg]", ylabel="Divergence speed [m/s]", xlims=[-60,0], legend=:top, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(ΛRange*180/π, UD, c=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!(ΛRange*180/π, UDRef, c=:black, ms=ms, label="Hodges and Pierce (2011)")
display(plt_UD)
savefig(string(absPath,"/sweptWingBendingDivergenceSweepRange_UD.pdf"))

# Divergence speed error vs. sweep angle
plt_eps = plot(xlabel="Sweep angle [deg]", ylabel="Error in divergence speed [%]", xlims=[-60,0], tickfont=font(ts), guidefont=font(fs))
plot!(ΛRange*180/π, (UD./UDRef .- 1)*100, c=:black, lw=lw, ls=:solid, label=false)
display(plt_eps)
savefig(string(absPath,"/sweptWingBendingDivergenceSweepRange_error.pdf"))

println("Finished sweptWingBendingDivergenceSweepRangePlotGenerator.jl")