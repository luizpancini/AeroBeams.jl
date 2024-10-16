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
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
msw = 0
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[URange[1],URange[end]])
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBflutter_AoA.pdf"))

# Trim thrust force vs airspeed
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
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-30,5],ylims=[0,120])
scatter!([NaN],[NaN], c=:black, shape=:circle, ms=ms, msw=msw, label="AeroBeams")
scatter!([NaN],[NaN], c=:black, shape=:utriangle, ms=ms, msw=msw, label="UM/NAST")
for mode in 1:nModes
    scatter!(dampsRef[mode+1,:], 2π*freqsRef[mode+1,:], c=:black, shape=:utriangle, ms=ms, msw=msw, label=false)
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
display(plt4)
savefig(string(absPath,"/BWBflutter_rootlocus.pdf"))

# V-g-f
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,120])
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw,  label=false)
end
plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.1])
for mode in 1:nModes
    scatter!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/BWBflutter_Vgf.pdf"))

println("Finished BWBflutterPlotGenerator.jl")