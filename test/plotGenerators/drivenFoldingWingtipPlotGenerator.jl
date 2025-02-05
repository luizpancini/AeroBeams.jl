using Plots

# Run the script
include("../examples/drivenFoldingWingtip.jl")

# Set paths
relPath = "/test/outputs/figures/drivenFoldingWingtip"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,view=(0,0),save=true,savePath=string(relPath,"/drivenFoldingWingtip_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# OOP position
plt_u3 = plot((r_n1.+u1)/x1_n[end], (r_n3.+u3)/x1_n[end], lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position", aspect_ratio=:equal)
display(plt_u3)
savefig(string(absPath,"/drivenFoldingWingtip_u3.pdf"))

# IP position
plt_u2 = plot((r_n1.+u1)/x1_n[end], (r_n2.+u2)/x1_n[end], lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized in-plane position", aspect_ratio=:equal)
display(plt_u2)
savefig(string(absPath,"/drivenFoldingWingtip_u2.pdf"))

# p1
plt3 = plot(x1_n/x1_n[end], p1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt3)
savefig(string(absPath,"/drivenFoldingWingtip_p1.pdf"))

# p2
plt2 = plot(x1_n/x1_n[end], p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/drivenFoldingWingtip_p2.pdf"))

# p3
plt2 = plot(x1_n/x1_n[end], p3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_3\$")
display(plt2)
savefig(string(absPath,"/drivenFoldingWingtip_p3.pdf"))

# Shear force
plt3 = plot(x1_n/x1_n[end], F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt3)
savefig(string(absPath,"/drivenFoldingWingtip_F3.pdf"))

# Bending moment
plt4 = plot(x1_n/x1_n[end], M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/drivenFoldingWingtip_M2.pdf"))

println("Finished drivenFoldingWingtipPlotGenerator.jl")