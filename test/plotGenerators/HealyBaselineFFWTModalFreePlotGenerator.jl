using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTModalFree.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTModalFree"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 6
msw = 0
gr()

# Mode shapes at highest airspeed
modesPlot = plot_mode_shapes(problem,scale=1,view=(15,15),frequencyLabel="frequency",legendPos=:outerright,save=true,savePath=string(relPath,"/HealyFFWTBeamModal_modeShapes.pdf"))
display(modesPlot)

println("Finished HealyBaselineFFWTModalFreePlotGenerator.jl")