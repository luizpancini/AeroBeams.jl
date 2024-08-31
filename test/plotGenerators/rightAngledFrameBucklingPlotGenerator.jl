using Plots

# Run the script
include("../examples/rightAngledFrameBuckling.jl")

# Set paths
relPath = "/test/outputs/figures/rightAngledFrameBuckling"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,view=(45,20),savePath=string(relPath,"/rightAngledFrameBuckling_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Plot normalized tip out-of-plane displacement over load steps
plt1 = plot(ylabel="\$F\$ [N]", xlabel="Tip \$u_2/L\$", title="Tip out-of-plane displacement")
plot!(tip_u2/L, σVector*F, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/rightAngledFrameBuckling_disp.pdf"))

println("Finished rightAngledFrameBucklingPlotGenerator.jl")