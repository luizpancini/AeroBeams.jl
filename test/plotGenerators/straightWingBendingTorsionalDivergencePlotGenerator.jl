using Plots, ColorSchemes

# Run the script
include("../examples/straightWingBendingTorsionalDivergence.jl")

# Plot configurations
lw = 2
ms = 4
gr()

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", ylims=[-50,5],legend=:topleft)
for i in eachindex(URange)
    scatter!(URange[i]*ones(length(dampingsNonOscillatory[i])), dampingsNonOscillatory[i], c=:black, ms=ms, msw=0, label=false)
end
display(plt_damp)
savefig(string(absPath,"/straightWingBendingTorsionalDivergence.pdf"))

println("Finished straightWingBendingTorsionalDivergence.jl")