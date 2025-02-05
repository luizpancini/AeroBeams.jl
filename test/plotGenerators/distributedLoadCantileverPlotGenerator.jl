using Plots

# Run the script
include("../examples/distributedLoadCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/distributedLoadCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,plotLimits=([-0.5,1],[-1,1],[-0.2,0.6]),save=true,savePath=string(relPath,"/distributedLoadCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
ts = 10
fs = 16
lw = 2
ms = 4
msw = 0
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="Tip generalized displacements", ylabel="\$q\$ [kip]",legend=:bottomright, tickfont=font(ts), guidefont=font(fs), legendfontsize=8)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-tip_u1/L, tip_u3/L, -tip_angle/π], colors)
    plot!(x, σVector*q, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1Ref[1,:], u3Ref[1,:], θRef[1,:]], [u1Ref[2,:], u3Ref[2,:], θRef[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/distributedLoadCantilever_disp.pdf"))

println("Finished distributedLoadCantileverPlotGenerator.jl")