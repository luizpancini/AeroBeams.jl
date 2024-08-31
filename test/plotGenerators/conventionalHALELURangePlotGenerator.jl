using Plots, ColorSchemes

# Run the script
include("../examples/conventionalHALELURange.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALELURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes of flexible aircraft at lowest airspeed
modesPlot = plot_mode_shapes(eigenProblem[1,1],nModes=5,scale=10,view=(30,30),legendPos=:outertop,save=true,savePath=string(relPath,"/conventionalHALELURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
lw = 2
ms = 5
msw = 0
lstyle = [:solid, :dash, :dot]
mshape = [:circle, :star, :utriangle]
gr()

# Root locus
plt0 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", ylims=[0,50])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=msw, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
display(plt0)
savefig(string(absPath,"/conventionalHALELURange_rootlocus.pdf"))

# Root locus
plt1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,0], ylims=[0,10])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=msw, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/conventionalHALELURange_rootlocus2.pdf"))

# Root locus (zoom)
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.1,0.1], ylims=[0,0.5])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=msw, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/conventionalHALELURange_rootlocuszoom.pdf"))

println("Finished conventionalHALELURangePlotGenerator.jl")