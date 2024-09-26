using Plots

# Run the script
include("../examples/PazyWingModal.jl")

# Build directory
mkpath(string(pwd(),"/test/outputs/figures/PazyWingModal"))

# Plot mode shapes
for i in [1,2]
    savePath = string("/test/outputs/figures/PazyWingModal/PazyWingModal_",i,".pdf")
    modesPlot = plot_mode_shapes(problem[i],scale=0.5,modalColorScheme=:rainbow,view=(30,30),frequencyLabel="frequency",save=true,savePath=savePath)
    display(modesPlot)
end

println("Finished PazyWingModalPlotGenerator.jl")