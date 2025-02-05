using Plots

# Run the script
include("../examples/spinupRobotArm.jl")

# Set paths
relPath = "/test/outputs/figures/spinupRobotArm"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="I",plotFrequency=5,showScale=false,timeStampPos=[0.1;-0.05;0],plotLimits=([-L,L],[-L,L],[0,L]),save=true,savePath=string(relPath,"/spinupRobotArm_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
ms = 5
msw = 0
gr()

# Normalized tip u₁
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_1/L\$")
plot!(t,u1_tip/L, c=:black, lw=lw, label="Flexible")
scatter!(t[1:5:end],u1_tip_rigid[1:5:end]/L, c=:blue, ms=ms, msw=msw, label="Rigid")
display(plt1)
savefig(string(absPath,"/spinupRobotArm_u1.pdf"))

# Normalized tip u₂
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_2/L\$")
plot!(t,u2_tip/L, c=:black, lw=lw, label="Flexible")
scatter!(t[1:5:end],u2_tip_rigid[1:5:end]/L, c=:blue, ms=ms, msw=msw, label="Rigid")
display(plt2)
savefig(string(absPath,"/spinupRobotArm_u2.pdf"))

# Root rotation
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root \$\\theta/\\pi\$")
plot!(t,θ3_root/π, c=:black, lw=lw, label="Numerical")
scatter!(t[1:20:end],θ(t[1:20:end])/π, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt3)
savefig(string(absPath,"/spinupRobotArm_theta.pdf"))

# Axial force
plot_time_outputs(problem,nodes=[(1,1)],elements=[1,nElements],nodalOutputs=["F1"],elementalOutputs=["F1"],save=true,saveFolder=string(relPath,"/"))

println("Finished spinupRobotArmPlotGenerator.jl")