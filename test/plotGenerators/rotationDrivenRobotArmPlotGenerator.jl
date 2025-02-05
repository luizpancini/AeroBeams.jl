using Plots

# Run the script
include("../examples/rotationDrivenRobotArm.jl")

# Set paths
relPath = "/test/outputs/figures/rotationDrivenRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=([-L,L],[0,L],[0,L]),save=true,savePath=string(relPath,"/rotationDrivenRobotArm_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,u1_tip/L, lw=lw, label="\$u_1/L\$")
plot!(t,u3_tip/L, lw=lw, label="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/rotationDrivenRobotArm_disp.pdf"))

println("Finished rotationDrivenRobotArmPlotGenerator.jl")