using Plots 

# Run the script
include("../examples/clampedFlaredHingedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/clampedFlaredHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(0,90),scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/clampedFlaredHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u2
plt1 = plot((r_n1.+u1)/x1[end], (r_n2.+u2)/x1[end], lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_2/L\$")
display(plt1)
savefig(string(absPath,"/clampedFlaredHingedBeam_u2.pdf"))

# u3
plt1 = plot((r_n1.+u1)/x1[end], (r_n3.+u3)/x1[end], lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedFlaredHingedBeam_u3.pdf"))

# p1
plt2 = plot(x1/L, p1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt2)
savefig(string(absPath,"/clampedFlaredHingedBeam_p1.pdf"))

# p2
plt3 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt3)
savefig(string(absPath,"/clampedFlaredHingedBeam_p2.pdf"))

# p3
plt4 = plot(x1/L, p3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_3\$")
display(plt4)
savefig(string(absPath,"/clampedFlaredHingedBeam_p3.pdf"))

# F3
plt5 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
display(plt5)
savefig(string(absPath,"/clampedFlaredHingedBeam_F3.pdf"))

# M2
plt6 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt6)
savefig(string(absPath,"/clampedFlaredHingedBeam_M2.pdf"))

println("Finished clampedFlaredHingedBeamPlotGenerator.jl")