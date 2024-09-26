using Plots

# Run the script
include("../examples/twoStoryFrame.jl")

# Set paths
relPath = "/test/outputs/figures/twoStoryFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=1,view=(45,30),legendPos=(0.3,0.1),frequencyLabel="frequency",save=true,savePath=string(relPath,"/twoStoryFrame_modeShapes.pdf"))
display(modesPlot)

println("Finished twoStoryFramePlotGenerator.jl")