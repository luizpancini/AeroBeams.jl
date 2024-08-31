using Plots

# Run the script
include("../examples/cantileverUnderSelfWeight.jl")

# Set paths
relPath = "/test/outputs/figures/cantileverUnderSelfWeight"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/cantileverUnderSelfWeight_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u1
plt0 = plot(x1/L, u1/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_1/L\$")
display(plt0)
savefig(string(absPath,"/cantileverUnderSelfWeight_u1.pdf"))

# u3
plt1 = plot(x1/L, u3/L, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(absPath,"/cantileverUnderSelfWeight_u3.pdf"))

# F3
plt2 = plot(x1/L, F3, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(absPath,"/cantileverUnderSelfWeight_F3.pdf"))

# M2
plt3 = plot(x1/L, M2, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(absPath,"/cantileverUnderSelfWeight_M2.pdf"))

println("Finished cantileverUnderSelfWeight.jl")