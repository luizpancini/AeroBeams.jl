using Plots

# Run the script
include("../examples/archUnderFollowerPressure.jl")

# Set paths
relPath = "/test/outputs/figures/archUnderFollowerPressure"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/archUnderFollowerPressure_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Plot normalized displacements over load steps
plt1 = plot(-mid_u3/R, σVector*λ, color=:black, lw=lw, xlabel="Midpoint \$-u_3/R\$", ylabel="\$\\lambda\$", label=false)
display(plt1)
savefig(string(absPath,"/archUnderFollowerPressure_disp.pdf"))

println("Finished archUnderFollowerPressurePlotGenerator.jl")