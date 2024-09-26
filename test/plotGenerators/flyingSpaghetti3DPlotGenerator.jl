using Plots

# Run the script
include("../examples/flyingSpaghetti3D.jl")

# Set paths
relPath = "/test/outputs/figures/flyingSpaghetti3D"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="I",plotFrequency=50,view=(30,30),plotLimits=[(0,2*L),(-L,L),(0,2*L)],save=true,savePath=string(relPath,"/flyingSpaghetti3D_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Nomalized tip displacements
labels = ["\$u_1/L\$" "\$u_2/L\$" "\$u_3/L\$"]
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u1_tip/L, u2_tip/L, u3_tip/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingSpaghetti3D_disp.pdf"))

# Nomalized tip angle
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$\\theta/(2\\pi)\$")
plot!(t,θ_tip/(2π), lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/flyingSpaghetti3D_angle.pdf"))

println("Finished flyingSpaghetti3DPlotGenerator.jl")