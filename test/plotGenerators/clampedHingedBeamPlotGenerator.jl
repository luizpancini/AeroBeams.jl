using Plots 

# Run the script
include("../examples/clampedHingedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/clampedHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1,showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/clampedHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot((r_n1.+u1)/x1[end], (r_n3.+u3)/x1[end], lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/clampedHingedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
display(plt3)
savefig(string(absPath,"/clampedHingedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt4)
savefig(string(absPath,"/clampedHingedBeam_M2.pdf"))

println("Finished clampedHingedBeamPlotGenerator.jl")