using Plots, ColorSchemes

# Run the script
include("../examples/cantileverWithTipTorsionalInertiaEigen.jl")

# Set paths
relPath = "/test/outputs/figures/cantileverWithTipTorsionalInertiaEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
gr()

# Plot mode shapes
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$p_1\$")
for m in 1:nModes
    plot!(x1/L, p1_modeShapes[m], lw=lw, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)
savefig(string(absPath,"/cantileverWithTipTorsionalInertiaEigen_p1.pdf"))

println("Finished cantileverWithTipTorsionalInertiaEigenPlotGenerator.jl")