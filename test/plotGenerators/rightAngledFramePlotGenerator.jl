using Plots

# Run the script
include("../examples/rightAngledFrame.jl")

# Set paths
relPath = "/test/outputs/figures/rightAngledFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,legendPos=:bottomleft,savePath=string(relPath,"/rightAngledFrame_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$-u_1/L, u_3/L, -\\theta/\\pi\$", ylabel="\$F\$ [N]", title="Tip generalized displacements")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-tip_u1/L, tip_u3/L, tip_angle/π], colors)
    plot!(x, σVector*F, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1Ref[1,:], u3Ref[1,:], θRef[1,:]], [u1Ref[2,:], u3Ref[2,:], θRef[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/rightAngledFrame_disp.pdf"))

println("Finished rightAngledFramePlotGenerator.jl")