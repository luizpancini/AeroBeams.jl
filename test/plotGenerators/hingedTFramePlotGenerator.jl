using Plots

# Run the script
include("../examples/hingedTFrame.jl")

# Set paths
relPath = "/test/outputs/figures/hingedTFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=300,view=(60,20),legendPos=(0.2,0.5),save=true,savePath=string(relPath,"/hingedTFrame_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
beamLabels = ["Beam 1" "Beam 2"]
forceLabels = ["\$F_3\$" "\$M_2\$"]
gr()

# u1
plt1 = plot([x1_beam1 x1_beam2], [u1_beam1 u1_beam2]/L, lw=lw, label=beamLabels, xlabel="\$x_1/L\$", ylabel="\$u_1/L\$")
display(plt1)
savefig(string(absPath,"/hingedTFrame_u1.pdf"))

# u3
plt2 = plot([x1_beam1 x1_beam2], [u3_beam1 u3_beam2]/L, lw=lw, label=beamLabels, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt2)
savefig(string(absPath,"/hingedTFrame_u3.pdf"))

# Internal loads
plt3 = plot(x1_beam1, [F3_beam1 M2_beam1], lw=lw, label=forceLabels, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N], \$M_2\$ [N.m]", title="Internal loads on beam 1")
display(plt3)
savefig(string(absPath,"/hingedTFrame_loads.pdf"))

println("Finished hingedTFramePlotGenerator.jl")