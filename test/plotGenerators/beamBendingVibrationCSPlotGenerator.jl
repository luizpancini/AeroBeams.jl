using Plots, ColorSchemes

# Run the script
include("../examples/beamBendingVibrationCS.jl")

# Set paths
relPath = "/test/outputs/figures/beamBendingVibrationCS"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
gr()

# Plot torsional mode shapes
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3\$")
for m in 1:nModes
    plot!(x1/L, u3_modeShapes[m], lw=lw, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)
savefig(string(absPath,"/beamBendingVibrationCS_u3.pdf"))

println("Finished beamBendingVibrationCSPlotGenerator.jl")