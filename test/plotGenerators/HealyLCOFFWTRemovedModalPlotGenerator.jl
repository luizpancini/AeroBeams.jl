using Plots

# Run the script
include("../examples/HealyLCOFFWTRemovedModal.jl")

# Build directory
mkpath(string(pwd(),"/test/outputs/figures/HealyLCOFFWTRemovedModal"))

# Plot all mode shapes
allModesPlot = plot_mode_shapes(problem,scale=1,lw=2,modeLabels=modeLabels,modalColorScheme=:rainbow,view=(45,15),plotAxes=false,save=true,savePath="/test/outputs/figures/HealyLCOFFWTRemovedModal/HealyLCOFFWTRemovedModal.pdf")
display(allModesPlot)

# Plot mode shapes one by one
modes2plot = 1:nModes
colors = cgrad(:rainbow, length(modes2plot), categorical=true)
for m in modes2plot
    p = plot_mode_shapes(problem,scale=1,lw=2,plotSteady=false,modes2plot=[m],modeLabels=[modeLabels[m]],modalColorScheme=palette([colors[m],colors[m]],2),view=(45,15),plotAxes=false,save=true,savePath=string("/test/outputs/figures/HealyLCOFFWTRemovedModal/HealyLCOFFWTRemovedModal_mode",m,".pdf"))
    display(p)
end

println("Finished HealyLCOFFWTRemovedModalPlotGenerator.jl")