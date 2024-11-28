using Plots 

# Run the script
include("../examples/biclampedHingedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/biclampedHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
scale = 1
deformationPlot = plot_steady_deformation(problem,scale=scale,showScale=true,scalePos=[0.3,0.1],save=true,savePath=string(relPath,"/biclampedHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/biclampedHingedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/biclampedHingedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt3)
savefig(string(absPath,"/biclampedHingedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/biclampedHingedBeam_M2.pdf"))

println("Finished biclampedHingedBeamPlotGenerator.jl")