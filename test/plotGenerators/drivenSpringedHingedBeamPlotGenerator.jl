using Plots

# Run the script
include("../examples/drivenSpringedHingedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/drivenSpringedHingedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,save=true,savePath=string(relPath,"/drivenSpringedHingedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# OOP position
plt_u3 = plot((r_n1.+u1)/L, (r_n3.+u3)/L, lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position", aspect_ratio=:equal)
display(plt_u3)
savefig(string(absPath,"/drivenSpringedHingedBeam_u3.pdf"))

# p2
plt2 = plot(x1_n/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/drivenSpringedHingedBeam_p2.pdf"))

# Shear force
plt_F3 = plot(x1_n/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
display(plt_F3)
savefig(string(absPath,"/drivenSpringedHingedBeam_F3.pdf"))

# Bending moment
plt_M2 = plot(x1_n/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt_M2)
savefig(string(absPath,"/drivenSpringedHingedBeam_M2.pdf"))

println("Finished drivenSpringedHingedBeamPlotGenerator.jl")