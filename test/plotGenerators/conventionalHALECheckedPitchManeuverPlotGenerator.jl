using Plots

# Run the script
include("../examples/conventionalHALECheckedPitchManeuver.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALECheckedPitchManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(dynamicProblem,refBasis="I",view=(30,15),plotDistLoads=false,plotFrequency=10,plotLimits=[(-20,20),(-10,250),(-20,20)],save=true,savePath=string(relPath,"/conventionalHALECheckedPitchManeuver_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Forward distance
plt1 = plot(xlabel="Time [s]", ylabel="Forward distance [m]")
plot!(t, Δu2, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALECheckedPitchManeuver_fdist.pdf"))

# Altitude
plt2 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALECheckedPitchManeuver_altitude.pdf"))

# AoA
plt3 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, wingAoA*180/π, c=:black, lw=lw, label="Wing")
plot!(t, htAoA*180/π, c=:blue, lw=lw, label="HT")
display(plt3)
savefig(string(absPath,"/conventionalHALECheckedPitchManeuver_wingAoA.pdf"))

# Airspeed
plt4 = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]")
plot!(t, airspeed, c=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/conventionalHALECheckedPitchManeuver_airspeed.pdf"))

println("Finished conventionalHALECheckedPitchManeuverPlotGenerator.jl")