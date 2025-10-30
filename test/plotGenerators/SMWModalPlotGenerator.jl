using Plots

# Run the script
include("../examples/SMWModal.jl")

# Set paths
relPath = "/test/outputs/figures/SMWModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=5,view=(30,30),legendPos=:best,save=true,savePath=string(relPath,"/SMWModal_modeShapes.pdf"))
display(modesPlot)

println("Finished SMWModalPlotGenerator.jl")