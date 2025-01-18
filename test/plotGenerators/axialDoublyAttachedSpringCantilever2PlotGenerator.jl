using Plots

# Run the script
include("../examples/axialDoublyAttachedSpringCantilever2.jl")

# Set paths
relPath = "/test/outputs/figures/axialDoublyAttachedSpringCantilever2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,save=true,savePath=string(relPath,"/axialDoublyAttachedSpringCantilever2_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Axial displacement
plt_u1 = plot(x1_n/L, u1/u₀, lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$u_1/u_0\$")
display(plt_u1)
savefig(string(absPath,"/axialDoublyAttachedSpringCantilever2_u1.pdf"))

# Axial force
plt_F1 = plot(x1_n/L, F1/F1[1], lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$F_1^{\\star}(x_1)/F_1^{\\star}(0)\$")
display(plt_F1)
savefig(string(absPath,"/axialDoublyAttachedSpringCantilever2_F1.pdf"))

println("Finished axialDoublyAttachedSpringCantilever2PlotGenerator.jl")