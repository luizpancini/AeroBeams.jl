using Plots

# Run the script
include("../examples/tipMomentCantileverDynamic.jl")

# Set paths
relPath = "/test/outputs/figures/tipMomentCantileverDynamic"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 20
fs = 24
lw = 2
fps = 5
pltFreq = round(Int,tf/Δt/20)
DPI = 300

# Animation
anim = plot_dynamic_deformation(problem,backendSymbol=:pyplot,plotAxes=false,plotFrequency=pltFreq,DPI=DPI,lw=lw,fps=fps,plotLimits=([-L/2,L],[0,L],[-3L/4,3L/4]),view=(0,0),save=true,savePath=string(relPath,"/tipMomentCantileverDynamic_deformation.gif"),displayProgress=true)
display(anim)

# Animated plot of curvature
gr()
tSpaced = t[1:pltFreq:end]
κ2Spaced = κ2[1:pltFreq:end]
plt = plot(xlabel="\$M_2\$", ylabel="\$\\kappa_2\$", title="\$\\kappa_2 = K_2 - k_2\$", xlims=[0,tf], ylims=[0, 2π], yticks=([0, 2π], ["\$0\$", "\$2π\$"]), xticks=false, tickfont=font(ts), guidefont=font(fs), titlefont=font(fs), dpi=DPI)
anim = @animate for (tt,timeNow) in enumerate(tSpaced)
    plot!(tSpaced[1:tt], κ2Spaced[1:tt], c=:black, label=false)
end
gif_handle = gif(anim, string(absPath,"/tipMomentCantileverDynamicPlotGenerator_kappa2.gif"), fps=fps)
display(gif_handle)

println("Finished tipMomentCantileverDynamicPlotGenerator.jl")