using Plots

# Run the script
include("../examples/doublePendulum.jl")

# Set paths
relPath = "/test/outputs/figures/doublePendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotFrequency=10,fps=60,scale=1,plotUndeformed=false,plotLimits=[(-L,L),(-L,0),(-L,L)],save=true,savePath=string(relPath,"/doublePendulum_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
labels = ["Hinge" "Tip"]
colors = [:black,:blue]
gr()

# Normalized tip u1 displacement
plt1 = plot(palette=colors, xlabel="\$t\$ [s]", ylabel="\$u_1/L\$ ")
plot!(t,[u1_hinge,u1_tip]/L, lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/doublePendulum_u1.pdf"))

# Normalized tip u3 displacement
plt2 = plot(palette=colors, xlabel="\$t\$ [s]", ylabel="\$u_3/L\$ ")
plot!(t,[u3_hinge,u3_tip]/L, lw=lw, label=labels)
display(plt2)
savefig(string(absPath,"/doublePendulum_u3.pdf"))

println("Finished doublePendulum.jl")