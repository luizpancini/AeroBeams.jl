using Plots

# Run the script
include("../examples/curvedCantileverDynamicFollower.jl")

# Set paths
relPath = "/test/outputs/figures/curvedCantileverDynamicFollower"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,scale=1,plotFrequency=5,plotLimits=([-L,L],[-L,L],[-L,L]),save=true,savePath=string(relPath,"/curvedCantileverDynamicFollower_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
dispLabels=["\$u_1/R\$" "\$u_2/R\$" "\$u_3/R\$"]
forceLabels=["\$F_1^*/F_0\$" "\$F_2^*/F_0\$" "\$F_3^*/F_0\$"]
momLabels=["\$M_1^*/(F_0R)\$" "\$M_2^*/(F_0R)\$" "\$M_3^*/(F_0R)\$"]
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip displacements")
plot!(t,[u1_tip/R,u2_tip/R,u3_tip/R], lw=lw, label=dispLabels)
display(plt1)
savefig(string(absPath,"/curvedCantileverDynamicFollower_disp.pdf"))

# Nomalized root forces
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Root forces")
plot!(t,[F1_root/F₀,F2_root/F₀,F3_root/F₀], lw=lw, label=forceLabels)
display(plt2)
savefig(string(absPath,"/curvedCantileverDynamicFollower_forces.pdf"))

# Nomalized root moments
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root moments")
plot!(t,[M1_root/(F₀*R),M2_root/(F₀*R),M3_root/(F₀*R)], lw=lw, label=momLabels)
display(plt3)
savefig(string(absPath,"/curvedCantileverDynamicFollower_moments.pdf"))

println("Finished curvedCantileverDynamicFollowerPlotGenerator.jl")