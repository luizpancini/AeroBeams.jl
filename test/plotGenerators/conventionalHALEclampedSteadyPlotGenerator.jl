using Plots

# Run the script
include("../examples/conventionalHALEclampedSteady.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALEclampedSteady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath=string(relPath,"/conventionalHALEclampedSteady_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan
L = maximum(x1_0)
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$ [% semispan]", xlims=[-1,1], ylims=[-10,30])
plot!(x1_def/L, x3_def/L*100, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALEclampedSteady_def.pdf"))

# Angle of attack over wingspan
plt2 = plot(xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]", xlims=[-1,1])
plot!(x1_e/L, α_of_x1*180/pi, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALEclampedSteady_aoa.pdf"))

println("Finished conventionalHALEclampedSteadyPlotGenerator.jl")