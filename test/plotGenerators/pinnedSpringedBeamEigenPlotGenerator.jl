using Plots

# Run the script
include("../examples/pinnedSpringedBeamEigen.jl")

# Set paths
relPath = "/test/outputs/figures/pinnedSpringedBeamEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=0.5,legendPos=:topleft,save=true,savePath=string(relPath,"/pinnedSpringedBeamEigen_modeShapes.pdf"))
display(modesPlot)

println("Finished pinnedSpringedBeamEigenPlotGenerator.jl")