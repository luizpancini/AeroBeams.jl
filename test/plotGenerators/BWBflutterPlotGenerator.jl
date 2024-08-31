using Plots, ColorSchemes

# Run the script
include("../examples/BWBflutter.jl")

# Set paths
relPath = "/test/outputs/figures/BWBflutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at highest airspeed
modesPlot = plot_mode_shapes(eigenProblem,scale=1,view=(30,30),legendPos=:outerright,save=true,savePath=string(relPath,"/BWBflutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[URange[1],URange[end]])
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBflutter_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/BWBflutter_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[URange[1],URange[end]])
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/BWBflutter_delta.pdf"))

# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-20,5],ylims=[0,120])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/BWBflutter_rootlocus.pdf"))

# V-g-f
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,120])
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0,  label=false)
end
plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=[URange[1],URange[end]], ylims=[-10,5])
for mode in 1:nModes
    scatter!(URange, modeDampings[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/BWBflutter_Vgf.pdf"))

println("Finished BWBflutterPlotGenerator.jl")