using Plots

# Run the script
include("../examples/cantileverWithTipSpring.jl")

# Set paths
relPath = "/test/outputs/figures/cantileverWithTipSpring"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,scale=1e2,showScale=true,scalePos=[0.25,-0.45],save=true,savePath=string(relPath,"/distributedLoadCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/cantileverWithTipSpring_u3.pdf"))

# F3
plt2 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/cantileverWithTipSpring_uF.pdf"))

# M2
plt3 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/cantileverWithTipSpring_M2.pdf"))

println("Finished cantileverWithTipSpring.jl")