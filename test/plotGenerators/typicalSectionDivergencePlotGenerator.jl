using Plots, ColorSchemes

# Run the script
include("../examples/typicalSectionDivergence.jl")

# Set paths
relPath = "/test/outputs/figures/typicalSectionDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 3
ms = 2
gr()

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", ylims=[-50,5],legend=:topleft)
for i in eachindex(URange)
    for j in eachindex(dampingsNonOscillatory[i])
        scatter!([URange[i]], [dampingsNonOscillatory[i][j]], c=:black, ms=ms, msw=0, label=false)
    end
end
display(plt_damp)
savefig(string(absPath,"/typicalSectionDivergence.pdf"))

println("Finished typicalSectionDivergencePlotGenerator.jl")