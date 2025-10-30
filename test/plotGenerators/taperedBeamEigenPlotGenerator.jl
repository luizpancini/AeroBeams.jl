using Plots

# Run the script
include("../examples/taperedBeamEigen.jl")

# Set paths
relPath = "/test/outputs/figures/taperedBeamEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=0.5,legendPos=:topleft,save=true,savePath=string(relPath,"/taperedBeamEigen_modeShapes.pdf"))
display(modesPlot)

println("Finished taperedBeamEigenPlotGenerators.jl")