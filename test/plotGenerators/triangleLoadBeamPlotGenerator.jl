using Plots

# Run the script
include("../examples/triangleLoadBeam.jl")

# Set paths
relPath = "/test/outputs/figures/triangleLoadBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e3,save=true,savePath=string(relPath,"/triangleLoadBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
gr()
lw = 2

# Displacement
gr()
plt1 = plot(x1/L, u3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3\$ [m]")
display(plt1)
savefig(string(absPath,"/triangleLoadBeam_u3.pdf"))
# Force
plt2 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/triangleLoadBeam_F3.pdf"))
# Moment
plt3 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/triangleLoadBeam_M2.pdf"))

println("Finished triangleLoadBeam.jl")