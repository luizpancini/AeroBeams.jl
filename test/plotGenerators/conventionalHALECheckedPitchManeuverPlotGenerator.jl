using Plots

# Run the script
include("../examples/conventionalHALECheckedPitchManeuver.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALECheckedPitchManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(dynamicProblem,refBasis="I",view=(90,0),plotDistLoads=false,plotFrequency=10,plotLimits=[(-20,20),(-10,100),(-50,50)],save=true,savePath=string(relPath,"/conventionalHALECheckedPitchManeuver_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALECheckedPitchManeuver_altitude.pdf"))

# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, rootAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALECheckedPitchManeuver_rootAoA.pdf"))

println("Finished conventionalHALECheckedPitchManeuverPlotGenerator.jl")