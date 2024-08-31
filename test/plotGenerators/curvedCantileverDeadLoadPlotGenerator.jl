using Plots

# Run the script
include("../examples/curvedCantileverDeadLoad.jl")

# Set paths
relPath = "/test/outputs/figures/curvedCantileverDeadLoad"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/curvedCantileverDeadLoad_deformation.pdf"))
display(deformationPlot)

println("Finished curvedCantileverDeadLoad.jl")