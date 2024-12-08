using Plots 

# Run the script
include("../examples/clampedHingedSpringedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/clampedHingedSpringedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(90+α*180/π,15),showScale=false,save=true,savePath=string(relPath,"/clampedHingedSpringedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u2
plt0 = plot((r_n1.+u1)/x1[end], (r_n2.+u2)/x1[end], lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_2/L\$")
display(plt0)
savefig(string(absPath,"/clampedHingedSpringedBeam_u2.pdf"))

# u3
plt1 = plot((r_n1.+u1)/x1[end], (r_n3.+u3)/x1[end], lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/clampedHingedSpringedBeam_u3.pdf"))

# p2_b
plt2 = plot(x1/L, p2_b, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2¨+\$")
display(plt2)
savefig(string(absPath,"/clampedHingedSpringedBeam_p2.pdf"))

# ps2
plt21 = plot(x1/L, ps2_b, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2¨+\$ (scaled)")
display(plt21)
savefig(string(absPath,"/clampedHingedSpringedBeam_ps2.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
display(plt3)
savefig(string(absPath,"/clampedHingedSpringedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt4)
savefig(string(absPath,"/clampedHingedSpringedBeam_M2.pdf"))

println("Finished clampedHingedSpringedBeamPlotGenerator.jl")