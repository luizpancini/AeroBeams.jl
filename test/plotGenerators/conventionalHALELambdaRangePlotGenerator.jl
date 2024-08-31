using Plots, ColorSchemes

# Run the script
include("../examples/conventionalHALELambdaRange.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALELambdaRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes of flexible aircraft
modesPlot = plot_mode_shapes(eigenProblem[end],nModes=5,scale=10,view=(30,30),legendPos=:outertop,save=true,savePath=string(relPath,"/conventionalHALELambdaRange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 4
xtick_positions = [1,2,5,10,20,50]
xtick_labels = ["1","2","5","10","20","50"]
ytick_positions = [0,10,20,50]
ytick_labels = ["0","10","20","50"]
gr()

# Trim root angle of attack vs stiffness factor
plt1 = plot(xlabel="Stiffness factor", ylabel="Trim root AoA [deg]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALELambdaRange_AoA.pdf"))

# Trim propeller force vs stiffness factor
plt2 = plot(xlabel="Stiffness factor", ylabel="Trim thrust [N]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALELambdaRange_thrust.pdf"))

# Trim elevator deflection vs stiffness factor
plt3 = plot(xlabel="Stiffness factor", ylabel="Trim elevator deflection [deg]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/conventionalHALELambdaRange_delta.pdf"))

# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,1], ylims=[0,300], yscale=:log10)
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/conventionalHALELambdaRange_rootlocus.pdf"))

# Frequency and damping evolution
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[0.9,maximum(λRange)], ylims=[0,50], xticks=(xtick_positions, xtick_labels), yticks=(ytick_positions, ytick_labels), xscale=:log10, yscale=:log10)
for mode in 1:nModes
    scatter!(λRange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0,  label=false)
end
plt52 = plot(xlabel="Stiffness factor", ylabel="Damping [1/s]", xlims=[0.9,maximum(λRange)], ylims=[-5,2], xscale=:log10, xticks=(xtick_positions, xtick_labels))
for mode in 1:nModes
    scatter!(λRange, modeDampings[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/conventionalHALELambdaRange_lambda-g-f.pdf"))

println("Finished conventionalHALELambdaRangePlotGenerator.jl")