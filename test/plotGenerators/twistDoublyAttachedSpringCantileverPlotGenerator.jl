using Plots

# Run the script
include("../examples/twistDoublyAttachedSpringCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/twistDoublyAttachedSpringCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,save=true,savePath=string(relPath,"/twistDoublyAttachedSpringCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Twist rotation
plt_p1 = plot(x1_n/L, p1, lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt_p1)
savefig(string(absPath,"/twistDoublyAttachedSpringCantilever_p1.pdf"))

# Twisting moment
plt_M1 = plot(x1_n/L, M1/M₀, lw=lw, label=false, xlims=[0,1], ylims=[0,1], xlabel="\$x_1/L\$", ylabel="\$M_1^{\\star}/M_0\$")
display(plt_M1)
savefig(string(absPath,"/twistDoublyAttachedSpringCantilever_M1.pdf"))

println("Finished twistDoublyAttachedSpringCantileverPlotGenerator.jl")