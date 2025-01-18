using Plots, ColorSchemes

# Run the script
include("../examples/PazyFFWTsteadyURangeAoARangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyURangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
gr()

# Coast angle vs airspeed for each AoA
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", ylims=[-135,45], yticks=-135:45:45)
for (i,θ) in enumerate(θRange)
    plot!(URange, -ϕHinge[i,:], lw=lw, c=colors[i], label="\$\\theta = $(θ*180/π) \\degree\$", legendfontsize=12, legend=:topleft)
end
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyURangeAoARangeCoast_deltaphihinge.pdf"))

println("Finished PazyFFWTsteadyURangeAoARangeCoastPlotGenerator.jl")