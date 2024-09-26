using Plots

# Run the script
include("../examples/tangentiallyForcedArch.jl")

# Set paths
relPath = "/test/outputs/figures/tangentiallyForcedArch"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/tangentiallyForcedArch_deformation.pdf"))
display(deformationPlot)

# Plot configurations
gr()
lw = 2
ms = 4
msw = 0
labels = ["\$-u_1/R\$" "\$u_3/R\$" "\$-\\theta/(\\pi/2)\$"]
XLabel = "\$-u_1/R, u_3/R, -\\theta/(\\pi/2)\$"
colors = [:blue,:orange,:green]

# Plot normalized displacements over load steps
plt1 = plot(xlabel=XLabel, ylabel="\$F\$ [kN]", title="Tip generalized displacements",legend=:bottomright)
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-tip_u1/R, tip_u3/R, -tip_angle/(π/2)], colors)
    plot!(x, σVector*abs(F)/(1e3), c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1Ref[1,:], u3Ref[1,:], θRef[1,:]], [u1Ref[2,:]/1e3, u3Ref[2,:]/1e3, θRef[2,:]/1e3], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/tangentiallyForcedArch_disp.pdf"))

println("Finished tangentiallyForcedArch.jl")