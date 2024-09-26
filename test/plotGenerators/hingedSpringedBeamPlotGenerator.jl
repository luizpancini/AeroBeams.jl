using Plots

# Run the script
include("../examples/hingedSpringedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/hingedSpringedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e2,showScale=true,scalePos=[0.3,0.5],save=true,savePath=string(relPath,"/hingedSpringedBeam_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/hingedSpringedBeam_u3.pdf"))

# p2
plt2 = plot(x1/L, p2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(absPath,"/hingedSpringedBeam_p2.pdf"))

# F3
plt3 = plot(x1/L, F3/1e3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [kN]")
display(plt3)
savefig(string(absPath,"/hingedSpringedBeam_F3.pdf"))

# M2
plt4 = plot(x1/L, M2/1e3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [kN.m]")
display(plt4)
savefig(string(absPath,"/hingedSpringedBeam_M2.pdf"))

println("Finished hingedSpringedBeamPlotGenerator.jl")