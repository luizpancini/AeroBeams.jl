using Plots

# Run the script
include("../examples/tipSineLoadedCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/tipSineLoadedCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,scale=1e3,plotFrequency=5,fps=30,plotLimits=([0,L],[-L,L],[-0.5,0.5]),save=true,savePath=string(relPath,"/tipSineLoadedCantilever_deformation.gif"),displayProgress=true)

# Plot configurations
gr()
lw = 1

# Tip u3
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [mm]")
plot!(t,u3_tip*1e3, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/tipSineLoadedCantilever_u3.pdf"))

# Root F3 
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_3^*\$ [N]")
plot!(t,F3_root, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/tipSineLoadedCantilever_F3.pdf"))

# Root M2 
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm]")
plot!(t,M2_root, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/tipSineLoadedCantilever_M2.pdf"))

println("Finished tipSineLoadedCantileverPlotGenerator.jl")