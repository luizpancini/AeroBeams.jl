using Plots, ColorSchemes

# Run the script
include("../examples/conventionalHALEURange.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALEURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at last airspeed
modesPlot = plot_mode_shapes(eigenProblem,nModes=5,scale=10,view=(30,30),legendPos=:outertop,save=true,savePath=string(relPath,"/conventionalHALEURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[20,35], ylims=[0,15])
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALEURange_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,35])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALEURange_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,35])
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/conventionalHALEURange_delta.pdf"))

# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-16.0,1.0], ylims=[0,50])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/conventionalHALEURange_rootlocus.pdf"))

# Root locus (zoom)
plt5 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-1,1], ylims=[0,10])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt5)
savefig(string(absPath,"/conventionalHALEURange_rootlocuszoom.pdf"))

println("Finished conventionalHALEURangePlotGenerator.jl")