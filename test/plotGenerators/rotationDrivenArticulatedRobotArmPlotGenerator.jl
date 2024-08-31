using Plots

# Run the script
include("../examples/rotationDrivenArticulatedRobotArm.jl")

# Set paths
relPath = "/test/outputs/figures/rotationDrivenArticulatedRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=[(-L,L),(-L/2,L),(0,L)],save=true,savePath=string(relPath,"/rotationDrivenArticulatedRobotArm_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
labels = ["Tip" "Hinge"]
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="\$u_1/L\$")
plot!(t,[u1_tip/L, u1_hinge/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/rotationDrivenArticulatedRobotArm_u1.pdf"))
plt2 = plot(xlabel="\$t\$ [s]", ylabel="\$u_3/L\$")
plot!(t,[u3_tip/L, u3_hinge/L], lw=lw, label=labels)
display(plt2)
savefig(string(absPath,"/rotationDrivenArticulatedRobotArm_u3.pdf"))

println("Finished rotationDrivenArticulatedRobotArmPlotGenerator.jl")