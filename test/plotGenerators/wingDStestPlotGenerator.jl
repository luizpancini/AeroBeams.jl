using Plots, ColorSchemes

# Run the script
include("../examples/wingDStest.jl")

# Set paths
relPath = "/test/outputs/figures/wingDStest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,view=(30,30),plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),save=true,savePath=string(relPath,"/wingDStest_deformation.gif"),displayProgress=true)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(elemRangePlot)))
labels = ["Root element", "Tip element"]
ts = 10
fs = 16
lw = 2
ms = 3
gr()

# Pitch angle
range2plot = 1:50:length(t)
θprescribed = (a₀-a₁.+θ.(t))*180/π
plt1 = plot(xlabel="Cycles", ylabel="Pitch angle [deg]", tickfont=font(ts), guidefont=font(fs), legend=:top, legendfontsize=10)
scatter!(t[range2plot]/τ, θprescribed[range2plot], color=:black, ms=ms, label="Prescribed at the root")
for (i,e) in enumerate(elemRangePlot)
    plot!(t/τ, α[i]*180/π, color=colors[i], lw=lw, label=labels[i])
end
display(plt1)
savefig(string(absPath,"/wingDStest_alpha.pdf"))

# cn vs time
plt2 = plot(xlabel="Time [s]", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t, cn[i], color=colors[i], lw=lw, label=labels[i])
end
display(plt2)
savefig(string(absPath,"/wingDStest_cnt.pdf"))

# cm vs time
plt3 = plot(xlabel="Time [s]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t, cm[i], color=colors[i], lw=lw, label=labels[i])
end
display(plt3)
savefig(string(absPath,"/wingDStest_cmt.pdf"))

# ct vs time
plt4 = plot(xlabel="Time [s]", ylabel="\$c_t\$", tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t, ct[i], color=colors[i], lw=lw, label=labels[i])
end
display(plt4)
savefig(string(absPath,"/wingDStest_ctt.pdf"))

# cl vs α
plt5 = plot(xlabel="Pitch angle [deg]", ylabel="\$c_l\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=10)
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Exp. - McAlister et al (1982)")
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i]*180/π, cl[i], color=colors[i], lw=lw, label=labels[i])
end
display(plt5)
savefig(string(absPath,"/wingDStest_cla.pdf"))

# cm vs α
plt6 = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legend=:bottomleft, legendfontsize=12)
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label=false)
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i]*180/π, cm[i], color=colors[i], lw=lw, label=false)
end
display(plt6)
savefig(string(absPath,"/wingDStest_cma.pdf"))

# cd vs α
plt7 = plot(xlabel="Pitch angle [deg]", ylabel="\$c_d\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=12)
scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label=false)
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i]*180/π, cdrag[i], color=colors[i], lw=lw, label=false)
end
display(plt7)
savefig(string(absPath,"/wingDStest_cda.pdf"))

println("Finished wingDStestPlotGenerator.jl")