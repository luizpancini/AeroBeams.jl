using Plots 

# Run the script
include("../examples/clampedHingedBeamRotated.jl")

# Set paths
relPath = "/test/outputs/figures/clampedHingedBeamRotated"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,15),showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/clampedHingedBeamRotated_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3 (seen from x2-plane)
plt1 = plot((r_n1.+u1)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedBeamRotated_u3_x2plane.pdf"))

# u3 (seen from x1-plane)
plt1 = plot((r_n2.+u2)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_2/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedBeamRotated_u3_x1plane.pdf"))

# p2_b (local basis)
plt2 = plot(x1/L, p2_b, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$p_2^+\$")
display(plt2)
savefig(string(absPath,"/clampedHingedBeamRotated_p2_b.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$F_3^{\\star}\$ [N]")
display(plt3)
savefig(string(absPath,"/clampedHingedBeamRotated_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt4)
savefig(string(absPath,"/clampedHingedBeamRotated_M2.pdf"))

println("Finished clampedHingedBeamRotatedPlotGenerator.jl")