using Plots 

# Run the script
include("../examples/sweptBackFFWTWing.jl")

# Set paths
relPath = "/test/outputs/figures/sweptBackFFWTWing"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath=string(relPath,"/sweptBackFFWTWing_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot((r_n1.+u1)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$",aspect_ratio=:equal)
display(plt1)
savefig(string(absPath,"/sweptBackFFWTWing_u3.pdf"))

# u2
plt1 = plot((r_n1.+u1)/L, (r_n2.+u2)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_2/L\$",aspect_ratio=:equal)
display(plt1)
savefig(string(absPath,"/sweptBackFFWTWing_u2.pdf"))

# p1
pltp1 = plot(x1/L, p1, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$p_1\$")
display(pltp1)
savefig(string(absPath,"/sweptBackFFWTWing_p1.pdf"))

# p2
pltp2 = plot(x1/L, p2, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$p_2\$")
display(pltp2)
savefig(string(absPath,"/sweptBackFFWTWing_p2.pdf"))

# p3
pltp3 = plot(x1/L, p3, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$p_3\$")
display(pltp3)
savefig(string(absPath,"/sweptBackFFWTWing_p3.pdf"))

# F3
plt3 = plot(x1/L, F3, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$F_3^{\\star}\$ [N]")
display(plt3)
savefig(string(absPath,"/sweptBackFFWTWing_F3.pdf"))

# M2
plt4 = plot(x1/L, M2, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$M_2^{\\star}\$ [N.m]")
display(plt4)
savefig(string(absPath,"/sweptBackFFWTWing_M2.pdf"))

println("Finished sweptBackFFWTWingPlotGenerator.jl")