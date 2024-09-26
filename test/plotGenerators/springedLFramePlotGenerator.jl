using Plots

# Run the script
include("../examples/springedLFrame.jl")

# Set paths
relPath = "/test/outputs/figures/springedLFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e2,save=true,savePath=string(relPath,"/springedLFrame_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Plots
plt1 = plot(x1/L2, u3_b/L2, lw=lw, label=false, xlabel="\$x_1/L_2\$", ylabel="\$u_3^{+}/L_2\$")
display(plt1)
savefig(string(absPath,"/springedLFrame_disp.pdf"))

println("Finished springedLFrame.jl")