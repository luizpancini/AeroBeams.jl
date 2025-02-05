using Plots

# Run the script
include("../examples/BWBcheckedPitchManeuver.jl")

# Set paths
relPath = "/test/outputs/figures/BWBcheckedPitchManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(dynamicProblem,refBasis="A",view=(30,15),plotDistLoads=false,plotFrequency=10,plotLimits=([-5,5],[-5,5],[-5,5]),save=true,savePath=string(relPath,"/BWBcheckedPitchManeuver_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBcheckedPitchManeuver_altitude.pdf"))

# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, rootAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/BWBcheckedPitchManeuver_rootAoA.pdf"))

println("Finished BWBcheckedPitchManeuverPlotGenerator.jl")