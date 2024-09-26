using Plots

# Run the script
include("../examples/conventionalHALECheckedRollManeuver.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALECheckedRollManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(dynamicProblem,refBasis="I",view=(30,15),plotDistLoads=false,plotFrequency=10,plotLimits=[(-200,20),(-10,200),(-50,50)],save=true,savePath=string(relPath,"/conventionalHALECheckedRollManeuver_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Lateral distance
plt1 = plot(xlabel="Time [s]", ylabel="Lateral distance [m]")
plot!(t, Δu1, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALECheckedRollManeuver_ldist.pdf"))

# Forward distance
plt1 = plot(xlabel="Time [s]", ylabel="Forward distance [m]")
plot!(t, Δu2, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALECheckedRollManeuver_fdist.pdf"))

# Altitude
plt2 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALECheckedRollManeuver_altitude.pdf"))

# AoA
plt3 = plot(xlabel="Time [s]", ylabel="Angle of attack [deg]")
plot!(t, wingAoA*180/π, c=:black, lw=lw, label="Wing (center)")
plot!(t, leftWingtipAoA*180/π, c=:green, lw=lw, label="Wing (left tip)")
plot!(t, rightWingtipAoA*180/π, c=:red, lw=lw, label="Wing (right tip.)")
plot!(t, htAoA*180/π, c=:blue, lw=lw, label="HT (center)")
plot!(t, vtAoA*180/π, c=:cyan, lw=lw, label="VT (root)")
display(plt3)
savefig(string(absPath,"/conventionalHALECheckedRollManeuver_AoA.pdf"))

# Airspeed
plt4 = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]")
plot!(t, airspeed, c=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/conventionalHALECheckedRollManeuver_airspeed.pdf"))

println("Finished conventionalHALECheckedRollManeuverPlotGenerator.jl")