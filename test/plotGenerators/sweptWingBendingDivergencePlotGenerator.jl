using Plots, ColorSchemes

# Run the script
include("../examples/sweptWingBendingDivergence.jl")

# Set paths
relPath = "/test/outputs/figures/sweptWingBendingDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 4
gr()

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", ylims=[-50,5],legend=:topleft)
for i in eachindex(URange)
    scatter!(URange[i]*ones(length(dampingsNonOscillatory[i])), dampingsNonOscillatory[i], c=:black, ms=2, msw=0, label=false)
end
display(plt_damp)
savefig(string(absPath,"/sweptWingBendingDivergence.pdf"))

println("Finished sweptWingBendingDivergencePlotGenerator.jl")