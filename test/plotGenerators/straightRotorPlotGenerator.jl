# Run the script
include("../examples/straightRotor.jl")

# Load experimental values from Epps & Chandra (1996)
ωRPMExp = [0,250,500,750]
expFreqs = [5.2374e+00   7.4239e+00   1.1018e+01   1.5820e+01;
            3.2999e+01   3.4582e+01   4.0188e+01   4.7001e+01;
            9.0439e+01   9.1112e+01   9.7733e+01   1.0558e+02;
            1.7237e+02   1.7532e+02   1.7967e+02   1.8874e+02]

# Set colormap and labels
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, 4))
labels = ["1B", "2B", "3B", "4B"]

# Plot mode shapes
relPath = "/test/outputs/figures/straightRotor"
absPath = string(pwd(),relPath)
mkpath(absPath)
modesPlot = plot_mode_shapes(problem,scale=5,view=(30,30),frequencyLabel="frequency",save=true,savePath=string(relPath,"/straightRotor_modeShapes.pdf"))
display(modesPlot)

# Plot frequency over angular velocity for several modes
gr()
lw = 2
ms = 5
msw = 0
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