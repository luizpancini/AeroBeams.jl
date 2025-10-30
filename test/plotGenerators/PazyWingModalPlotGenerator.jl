using Plots

# Run the script
include("../examples/PazyWingModal.jl")

# Build directory
mkpath(string(pwd(),"/test/outputs/figures/PazyWingModal"))

# Plot static mode shapes
for i in [1,2]
    legendPos = i==1 ? (0.3,0.55) : :best
    savePath = string("/test/outputs/figures/PazyWingModal/PazyWingModal_",i,".pdf")
    modesPlot = plot_mode_shapes(problem[i],scale=1,modeLabels=["OOP-1","OOP-2","T-1","OOP-3","T-IP-1"],modalColorScheme=:rainbow,view=(45,15),plotAxes=false,legendPos=legendPos,save=true,savePath=savePath)
    display(modesPlot)
end

# # Plot mode shapes animation
# for i in [1,2]
#     legendPos = i==1 ? (0.3,0.6) : :best
#     savePath = string("/test/outputs/figures/PazyWingModal/PazyWingModal_",i,".gif")
#     modesAnim = plot_mode_shapes_animation(problem[i],matchModeFrequency=true,nFramesPerCycle=101,modeLabels=["OOP-1","OOP-2","T-1","OOP-3","T-IP-1"],modalColorScheme=:rainbow,view=(45,15),scale=1,plotAxes=false,legendPos=legendPos,legendFontSize=10,save=true,savePath=savePath,displayProgress=true)
#     display(modesAnim)
# end

println("Finished PazyWingModalPlotGenerator.jl")