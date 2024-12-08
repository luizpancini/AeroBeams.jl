using Plots

# Run the script
include("../examples/PazyFFWTsteadyCoast.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,view=(60,30),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteadyCoast_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# OOP displacement
plt1 = plot((r_n1.+u1)/x1_n[end], (r_n3.+u3)/x1_n[end], lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position")
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyCoast_u3.pdf"))

# p1
plt3 = plot(x1_n/x1_n[end], p1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt3)
savefig(string(absPath,"/PazyFFWTsteadyCoast_p1.pdf"))

# p2
plt2 = plot(x1_n/x1_n[end], p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyCoast_p2.pdf"))

# p3
plt2 = plot(x1_n/x1_n[end], p3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_3\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyCoast_p3.pdf"))

# Bending moment
plt4 = plot(x1_n/x1_n[end], M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt4)
savefig(string(absPath,"/PazyFFWTsteadyCoast_M2.pdf"))

# Angle of attack
plt5 = plot(x1_e/x1_n[end], α*180/π, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]")
display(plt5)
savefig(string(absPath,"/PazyFFWTsteadyCoast_aoa.pdf"))

# cn
plt6 = plot(x1_e/x1_n[end], cn, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$c_n\$ ")
display(plt6)
savefig(string(absPath,"/PazyFFWTsteadyCoast_cn.pdf"))

println("Finished PazyFFWTsteadyCoastPlotGenerator.jl")