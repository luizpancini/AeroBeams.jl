using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyModal.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot swept configuration
plt_swept_modeShapes = plot_mode_shapes(problem[end],view=(30,30),interactive=false,frequencyLabel="frequency",save=true,savePath=string(relPath,"/PazySweepAngleRangeModal_modeShapes.pdf"))
display(plt_swept_modeShapes)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Sweep angle [deg]", ylabel="Frequency [Hz]", ylims=[0,50])
for mode in 1:nModes
    plot!(ΛRange*180/π, modeFrequencies[mode], c=colors[mode], lw=lw,  label=false)
end
display(plt1)
savefig(string(absPath,"/sweptPazyModal.pdf"))

println("Finished sweptPazyModal.jl")