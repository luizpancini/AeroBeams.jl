using Plots

# Run the script
include("../examples/compositeCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/compositeCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1,plotFrequency=10,plotLimits=([0,L],[-5,5],[-5,5]),save=true,savePath=string(relPath,"/compositeCantilever_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip u1
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_1\$ [m] ")
plot!(t,u1_tip, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/compositeCantilever_u1.pdf"))

# Tip u2 
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_2\$ [m] ")
plot!(t,u2_tip, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/compositeCantilever_u2.pdf"))

# Tip u3
plt3 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [m] ")
plot!(t,u3_tip, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/compositeCantilever_u3.pdf"))

# Tip p1 
plt4 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$p_1\$")
plot!(t,p1_tip, c=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/compositeCantilever_p1.pdf"))

# Tip p2 
plt5 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$p_2\$")
plot!(t,p2_tip, c=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/compositeCantilever_p2.pdf"))

# Tip p3
plt6 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$p_3\$")
plot!(t,p3_tip, c=:black, lw=lw, label=false)
display(plt6)
savefig(string(absPath,"/compositeCantilever_p3.pdf"))

# Root F1 
plt7 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_1^*\$ [N]")
plot!(t,F1_root, c=:black, lw=lw, label=false)
display(plt7)
savefig(string(absPath,"/compositeCantilever_F1.pdf"))

# Root F2
plt8 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_2^*\$ [N]")
plot!(t,F2_root, c=:black, lw=lw, label=false)
display(plt8)
savefig(string(absPath,"/compositeCantilever_F2.pdf"))

# Root F3 
plt9 = plot(xlabel="\$t\$ [s]", ylabel="Root \$F_3^*\$ [N]")
plot!(t,F3_root, c=:black, lw=lw, label=false)
display(plt9)
savefig(string(absPath,"/compositeCantilever_F3.pdf"))

# Root M1 
plt10 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_1^*\$ [Nm]")
plot!(t,M1_root, c=:black, lw=lw, label=false)
display(plt10)
savefig(string(absPath,"/compositeCantilever_M1.pdf"))

# Root M2 
plt11 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm]")
plot!(t,M2_root, c=:black, lw=lw, label=false)
display(plt11)
savefig(string(absPath,"/compositeCantilever_M2.pdf"))

# Root M3 
plt12 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_3^*\$ [Nm]")
plot!(t,M3_root, c=:black, lw=lw, label=false)
display(plt12)
savefig(string(absPath,"/compositeCantilever_M3.pdf"))

println("Finished compositeCantileverPlotGenerator.jl")