using Plots, ColorSchemes

# Run the script
include("../examples/wingDStest.jl")

# Set paths
relPath = "/test/outputs/figures/wingDStest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,view=(30,30),plotLimits=[(0,L),(-L/2,L/2),(-L/2,L/2)],save=true,savePath=string(relPath,"/wingDStest_deformation.gif"),displayProgress=true)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(elemRangePlot)))
labels = string.(elemRangePlot)
lw = 2
ms = 3
gr()

# Range of last cycle
indBeginLastCycle = argmin(abs.(t .- τ))
rangeLastCycle = indBeginLastCycle:length(t)

# Pitch angle
plt1 = plot(xlabel="Time [s]", ylabel="Pitch angle [deg]")
for (i,e) in enumerate(elemRangePlot)
    plot!(t, α[i]*180/π, color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt1)
savefig(string(absPath,"/wingDStest_alpha.pdf"))

# cn vs time
plt2 = plot(xlabel="Time [s]", ylabel="\$c_n\$")
for (i,e) in enumerate(elemRangePlot)
    plot!(t, cn[i], color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt2)
savefig(string(absPath,"/wingDStest_cnt.pdf"))

# cm vs time
plt3 = plot(xlabel="Time [s]", ylabel="\$c_m\$")
for (i,e) in enumerate(elemRangePlot)
    plot!(t, cm[i], color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt3)
savefig(string(absPath,"/wingDStest_cmt.pdf"))

# ct vs time
plt4 = plot(xlabel="Time [s]", ylabel="\$c_t\$")
for (i,e) in enumerate(elemRangePlot)
    plot!(t, ct[i], color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt4)
savefig(string(absPath,"/wingDStest_ctt.pdf"))

# cl vs α
plt5 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
plot!([NaN], [NaN], color=:black, lw=lw, label="AeroBeams")
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i][rangeLastCycle]*180/π, cl[i][rangeLastCycle], color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt5)
savefig(string(absPath,"/wingDStest_cla.pdf"))

# cm vs α
plt6 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$", legend=:bottomleft)
plot!([NaN], [NaN], color=:black, lw=lw, label="AeroBeams")
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i][rangeLastCycle]*180/π, cm[i][rangeLastCycle], color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt6)
savefig(string(absPath,"/wingDStest_cma.pdf"))

# cd vs α
plt7 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
plot!([NaN], [NaN], color=:black, lw=lw, label="AeroBeams")
scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i][rangeLastCycle]*180/π, cdrag[i][rangeLastCycle], color=colors[i], lw=lw, label=string("Element ",labels[i]))
end
display(plt7)
savefig(string(absPath,"/wingDStest_cda.pdf"))

println("Finished wingDStestPlotGenerator.jl")