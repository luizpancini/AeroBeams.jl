using Plots

# Run the script
include("../examples/PazyWingModal.jl")

# Build directory
mkpath(string(pwd(),"/test/outputs/figures/PazyWingModal"))

# Plot mode shapes
for i in [1,2]
    legendPos = i==1 ? (0.3,0.55) : :best
    savePath = string("/test/outputs/figures/PazyWingModal/PazyWingModal_",i,".pdf")
    modesPlot = plot_mode_shapes(problem[i],scale=1,lw=2,modalColorScheme=:rainbow,view=(45,15),plotAxes=false,legendPos=legendPos,frequencyLabel="frequency",save=true,savePath=savePath)
    display(modesPlot)
end

println("Finished PazyWingModalPlotGenerator.jl")