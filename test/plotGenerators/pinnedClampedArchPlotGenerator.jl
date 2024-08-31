using Plots

# Run the script
include("../examples/pinnedClampedArch.jl")

# Set paths
relPath = "/test/outputs/figures/pinnedClampedArch"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/pinnedClampedArch_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
x = [-u1_atForce/R, -u3_atForce/R]
labels = ["\$-u_1/R\$" "\$-u_3/R\$"]
colors = [:blue,:orange]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$-u_1/L, -u_3/L,\$", ylabel="\$\\lambda\$", title="Displacements at point of force application",legend=:bottomright)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="DaDeppo & Schmidt (1975)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*λ, lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]], [u1Ref[2,:],u3Ref[2,:]], palette=colors,ms=ms, msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/pinnedClampedArch_disp.pdf"))

println("Finished pinnedClampedArchPlotGenerator.jl")