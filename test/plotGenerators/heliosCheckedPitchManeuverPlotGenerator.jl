using Plots

# Run the script
include("../examples/heliosCheckedPitchManeuver.jl")

# Set paths
relPath = "/test/outputs/figures/heliosCheckedPitchManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,refBasis="I",view=(30,15),plotDistLoads=false,plotFrequency=5,plotLimits=([-30,30],[0,600],[-20,20]),save=true,savePath=string(relPath,"/heliosCheckedPitchManeuver_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/heliosCheckedPitchManeuver_altitude.pdf"))

# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, rootAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/heliosCheckedPitchManeuver_rootAoA.pdf"))

println("Finished heliosCheckedPitchManeuverPlotGenerator.jl")