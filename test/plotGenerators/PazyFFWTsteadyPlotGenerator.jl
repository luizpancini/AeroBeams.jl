using Plots

# Run the script
include("../examples/PazyFFWTsteady.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteady_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# OOP displacement
plt1 = plot((r_n1.+u1_of_x1)/x1[end], (r_n3.+u3_of_x1)/x1[end], aspect_ratio=:equal, lw=lw, label=false, xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position")
display(plt1)
savefig(string(absPath,"/PazyFFWTsteady_u3.pdf"))

# Bending angle
plt2 = plot(x1/x1[end], p2_of_x1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/PazyFFWTsteady_p2.pdf"))

# Bending moment
plt3 = plot(x1/x1[end], M2_of_x1, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/PazyFFWTsteady_M2.pdf"))

println("Finished PazyFFWTsteadyPlotGenerator.jl")