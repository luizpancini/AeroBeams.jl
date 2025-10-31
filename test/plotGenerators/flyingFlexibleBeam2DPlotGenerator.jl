using Plots

# Run the script
include("../examples/flyingFlexibleBeam2D.jl")

# Set paths
relPath = "/test/outputs/figures/flyingFlexibleBeam2D"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,refBasis="I",plotFrequency=10,plotLimits=([0,3*L],[-L,0],[-L/2,L/2]),save=true,savePath=string(relPath,"/flyingFlexibleBeam2D_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
labels = ["\$u_1/L\$" "\$u_3/L\$"]
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u1_tip/L, u3_tip/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingFlexibleBeam2D_disp.pdf"))

println("Finished flyingFlexibleBeam2DPlotGenerator.jl")