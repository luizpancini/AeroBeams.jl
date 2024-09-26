using Plots

# Run the script
include("../examples/flyingSpaghetti2D.jl")

# Set paths
relPath = "/test/outputs/figures/flyingSpaghetti2D"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="I",plotFrequency=10,plotLimits=[(0,2*L),(-L,0),(-L/2,L/2)],save=true,savePath=string(relPath,"/flyingSpaghetti2D_deformation.gif"),displayProgress=true)

# Plot configurations
labels = ["\$u_1/L\$" "\$u_3/L\$"]
lw = 2
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u1_tip/L, u3_tip/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingSpaghetti2D_disp.pdf"))

println("Finished flyingSpaghetti2DPlotGenerator.jl")