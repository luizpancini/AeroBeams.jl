using Plots

# Run the script
include("../examples/axialTractionCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/axialTractionCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,scale=1e4,plotLimits=[(0,2*L),(-0.1,0.1),(0,1)],save=true,savePath=string(relPath,"/axialTractionCantilever_deformation.gif"))

# Plot configurations
lw = 2
ms = 5
colors = [:blue,:orange]
labels = ["\$t\$ = 0.8 s" "\$t\$ = 1.0 s"]
gr()

# Axial displacement at time instances 0.8 s and 1.0 s
plt1 = plot(xlabel="\$x_1\$ [m]", ylabel="\$u_1\$ [mm]")
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Reddy (2005)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x1,u1_08*1e3, c=colors[1], lw=lw, label=false)
plot!(x1,u1_10*1e3, c=colors[2], lw=lw, label=false)
scatter!(x1_ref,u1_08_ref, c=colors[1], ms=ms, msw=msw, label=false)
scatter!(x1_ref,u1_10_ref, c=colors[2], ms=ms, msw=msw, label=false)
display(plt1)
savefig(string(absPath,"/axialTractionCantilever_u1.pdf"))

println("Finished axialTractionCantileverPlotGenerator.jl")