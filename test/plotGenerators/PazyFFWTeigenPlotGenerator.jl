using Plots

# Run the script
include("../examples/PazyFFWTeigen.jl")

# Set paths
mkpath(string(pwd(),"/test/outputs/figures/PazyFFWTeigen"))

# Mode shapes
modesPlot = plot_mode_shapes(problem,scale=0.1,view=(30,30),frequencyLabel="frequency",save=true,savePath="/test/outputs/figures/PazyFFWTeigen/PazyFFWTeigen_modeShapes.pdf")
display(modesPlot)

println("Finished PazyFFWTeigenPlotGenerator.jl")