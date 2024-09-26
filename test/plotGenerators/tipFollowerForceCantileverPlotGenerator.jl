using Plots, DelimitedFiles

# Run the script
include("../examples/tipFollowerForceCantilever.jl")

# Set paths
relPath = "/test/outputs/figures/tipFollowerForceCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomright,save=true,savePath=string(relPath,"/tipFollowerForceCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = [:blue,:green]
labels = ["\$-u_1/L\$" "\$u_3/L\$"]
lw = 2
ms = 3
msw = 0
gr()

# Plot normalized tip displacements over load steps
plt1 = plot(xlabel="\$F\$ [kip]", ylabel="\$-u_1/L, u_3/L\$", title="Tip generalized displacements", xticks=collect(0:10:F), yticks=collect(-0.6:0.2:1.2))
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1986)")
for i=1:2
    plot!([NaN], [NaN], lc=colors[i], m=colors[i],  lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(σVector*F/(1e3), [-tip_u1/L, tip_u3/L], palette=colors,  lw=lw, label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]], [u1Ref[2,:]/L,u3Ref[2,:]/L], palette=colors, ms=ms, msw=msw, label=false)
display(plt1)
savefig(string(absPath,"/tipFollowerForceCantilever_disp.pdf"))

println("Finished tipFollowerForceCantileverPlotGenerator.jl")