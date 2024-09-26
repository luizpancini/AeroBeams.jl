using Plots 

# Run the script
include("../examples/momentDrivenRobotArm.jl")

# Set paths
relPath = "/test/outputs/figures/momentDrivenRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotFrequency=1,showScale=false,timeStampPos=[0.15;-0.05;0],plotLimits=[(-L,L),(-L,L),(0,L)],save=true,savePath=string(relPath,"/momentDrivenRobotArm_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,u1_tip/L, lw=lw, label="\$u_1/L\$")
plot!(t,u3_tip/L, lw=lw, label="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/momentDrivenRobotArm_disp.pdf"))

println("Finished momentDrivenRobotArmPlotGenerator.jl")