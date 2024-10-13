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
lw = 2
ms = 3
gr()

# Range of last cycle
indBeginLastCycle = argmin(abs.(t .- τ))
rangeLastCycle = indBeginLastCycle:length(t)

# Pitch angle
plt1 = plot(xlabel="Time [s]", ylabel="Pitch angle [deg]")
plot!(t, α*180/π, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/wingDStest_alpha.pdf"))

# Normal relative wind acceleration
plt2 = plot(xlabel="Time [s]", ylabel="\$\\dot{V}_3\$ [m/s^2]")
plot!(t, Vdot3, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/wingDStest_Vdot3.pdf"))

# cn vs time
plt3 = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/wingDStest_cnt.pdf"))

# cm vs time
plt4 = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/wingDStest_cmt.pdf"))

# ct vs time
plt5 = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/wingDStest_ctt.pdf"))

# cl vs α
plt6 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
plot!(α[rangeLastCycle]*180/π, cl[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt6)
savefig(string(absPath,"/wingDStest_cla.pdf"))

# cm vs α
plt7 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$", legend=:bottomleft)
plot!(α[rangeLastCycle]*180/π, cm[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt7)
savefig(string(absPath,"/wingDStest_cma.pdf"))

# cd vs α
plt8 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
plot!(α[rangeLastCycle]*180/π, cdrag[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt8)
savefig(string(absPath,"/wingDStest_cda.pdf"))

println("Finished wingDStestPlotGenerator.jl")