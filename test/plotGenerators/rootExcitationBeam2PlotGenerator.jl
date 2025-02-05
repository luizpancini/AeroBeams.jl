using Plots

# Run the script
include("../examples/rootExcitationBeam2.jl")

# Set paths
relPath = "/test/outputs/figures/rootExcitationBeam2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotFrequency=5,scalePos=[1;-0.05;0],timeStampPos=[1;-0.15;0],plotLimits=([-0.1,0.1],[0,L],[0,L]),plotDistLoads=false,save=true,savePath=string(relPath,"/rootExcitationBeam2_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
ms = 3
gr()

# Normalized V3/V over the beam, over time
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$V_3/V_b\$")
for i in 1:length(t)
    plot!(x1/L,V3[i]/V, c=:black, lw=lw, label=false)
end
display(plt1)
savefig(string(absPath,"/rootExcitationBeam2_V3oV.pdf"))

# Normalized u3 (basis b) at the root and at the tip
plt2 = plot()
plot!(t/T,u3b_root/A, c=:auto, lw=lw, label="Root", xlabel="\$t/T\$")
plot!(t/T,u3b_tip/A, c=:auto, lw=lw, label="Tip", xlabel="\$t/T\$", ylabel="\$u_3^{+}/A\$")
display(plt2)
savefig(string(absPath,"/rootExcitationBeam2_u3.pdf"))

# Normalized V3 at the root and at the tip
plt3 = plot()
plot!(t/T,V3_root/V, c=:auto, lw=lw, label="Root", xlabel="\$t/T\$")
plot!(t/T,V3_tip/V, c=:auto, lw=lw, label="Tip", xlabel="\$t/T\$", ylabel="\$V_3^{*}/A\$")
display(plt3)
savefig(string(absPath,"/rootExcitationBeam2_V3.pdf"))

println("Finished rootExcitationBeam2PlotGenerator.jl")