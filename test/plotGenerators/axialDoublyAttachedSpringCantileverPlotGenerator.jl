using Plots

# Run the script
include("../examples/axialDoublyAttachedSpringCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/axialDoublyAttachedSpringCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,save=true,savePath=string(relPath,"/axialDoublyAttachedSpringCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Axial displacement
plt_u1 = plot(x1_n/L, u1/L, lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$u_1/L\$")
display(plt_u1)
savefig(string(absPath,"/axialDoublyAttachedSpringCantilever_u1.pdf"))

# Axial force
plt_F1 = plot(x1_n/L, F1/F₀, lw=lw, label=false, xlims=[0,1], ylims=[0,1], xlabel="\$x_1/L\$", ylabel="\$F_1^{\\star}/F_0\$")
display(plt_F1)
savefig(string(absPath,"/axialDoublyAttachedSpringCantilever_F1.pdf"))

println("Finished axialDoublyAttachedSpringCantileverPlotGenerator.jl")