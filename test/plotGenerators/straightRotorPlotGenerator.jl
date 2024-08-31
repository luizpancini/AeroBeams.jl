using Plots

# Run the script
include("../examples/straightRotor.jl")

# Set paths
relPath = "/test/outputs/figures/straightRotor"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=5,view=(30,30),frequencyLabel="frequency",save=true,savePath=string(relPath,"/straightRotor_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
gr()
lw = 2
ms = 5
msw = 0
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, 4))
labels = ["1B", "2B", "3B", "4B"]

# Plot frequency over angular velocity for several modes
plt1 = plot(xlabel="Angular velocity [rpm]", ylabel="Frequency [Hz]", title="Bending modes", xticks=ωRPMExp, yticks=collect(0:20:200))
modes = [1,2,5,6]
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
for (m,mode) in enumerate(modes)  
    numFreqsMode = [numFreqs[j][mode]/(2*π) for j in 1:length(ωRPM)]
    plot!(ωRPM,numFreqsMode, lc=colors[m], lw=lw, label=false)
    scatter!(ωRPMExp,expFreqs[m,:], mc=colors[m], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[m], m=colors[m], lw=lw, ms=ms, msw=msw, label=labels[m])
end
plot!(legend=(0.15,0.8))
display(plt1)
savefig(string(absPath,"/straightRotor_freqVsOmega.pdf"))

println("Finished straightRotorPlotGenerator.jl")