using Plots, ColorSchemes

# Run the script
include("../examples/sweptWingTorsionalDivergence.jl")

# Set paths
relPath = "/test/outputs/figures/sweptWingTorsionalDivergence"
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

println("Finished sweptWingTorsionalDivergencePlotGenerator.jl")