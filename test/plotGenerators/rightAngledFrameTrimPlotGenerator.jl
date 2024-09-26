using Plots

# Run the script
include("../examples/rightAngledFrameTrim.jl")

# Set paths
relPath = "/test/outputs/figures/rightAngledFrameTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomleft,save=true,savePath=string(relPath,"/rightAngledFrameTrim_deformation.pdf"))
display(deformationPlot)