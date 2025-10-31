using Plots

# Run the script
include("../examples/HealyLCOFFWTLockedModal.jl")

# Build directory
mkpath(string(pwd(),"/test/outputs/figures/HealyLCOFFWTLockedModal"))

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=1,lw=2,modeLabels=modeLabels,modalColorScheme=:rainbow,view=(45,15),plotAxes=false,save=true,savePath="/test/outputs/figures/HealyLCOFFWTLockedModal/HealyLCOFFWTLockedModal.pdf")
display(modesPlot)

println("Finished HealyLCOFFWTLockedModalPlotGenerator.jl")