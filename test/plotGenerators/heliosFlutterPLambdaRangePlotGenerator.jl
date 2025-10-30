using Plots, ColorSchemes

# Run the script
include("../examples/heliosFlutterPLambdaRange.jl")

# Set paths
relPath = "/test/outputs/figures/heliosFlutterPLambdaRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes of flexible aircraft at highest payload
modesPlot = plot_mode_shapes(eigenProblem[1,end],element2centralize=15,scale=10,view=(30,30),legendPos=:outerright,nModes=6,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
lw = 2
ms = 3
lstyle = [:solid, :dash, :dot]
mshape = [:circle, :star, :utriangle]
gr()

# Root locus 
plt1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/heliosFlutterPLambdaRange_rootlocus.pdf"))

# Root locus (phugoid zoom)
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.1,0.2], ylims=[0,0.6])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/heliosFlutterPLambdaRange_Phzoom.pdf"))

# Root locus (dutch roll zoom)
plt3 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.6,0], ylims=[0,0.6])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/heliosFlutterPLambdaRange_DRzoom.pdf"))

# Root locus (short-period zoom)
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,-2.5], ylims=[0,5])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt4)
savefig(string(absPath,"/heliosFlutterPLambdaRange_SPzoom.pdf"))

println("Finished heliosFlutterPLambdaRangePlotGenerator.jl")