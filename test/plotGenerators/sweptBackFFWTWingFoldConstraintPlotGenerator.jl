using Plots 

# Run the script
include("../examples/sweptBackFFWTWingFoldConstraint.jl")

# Set paths
relPath = "/test/outputs/figures/sweptBackFFWTWingFoldConstraint"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(Λ*180/pi,30),showScale=false,save=true,savePath=string(relPath,"/sweptBackFFWTWingFoldConstraint_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot((r_n1.+u1)/L, (r_n3.+u3)/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
display(plt1)
savefig(string(absPath,"/sweptBackFFWTWingFoldConstraint_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="Normalized arclength", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/sweptBackFFWTWingFoldConstraint_p2.pdf"))

println("Finished sweptBackFFWTWingFoldConstraintPlotGenerator.jl")