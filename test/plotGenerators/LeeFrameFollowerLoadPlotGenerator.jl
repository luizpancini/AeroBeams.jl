using Plots

# Run the script
include("../examples/LeeFrameFollowerLoad.jl")

# Set paths
relPath = "/test/outputs/figures/LeeFrameFollowerLoad"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomright,save=true,savePath=string(relPath,"/LeeFrameFollowerLoad_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 3
x = [u1_atForce/L, -u3_atForce/L]
labels = ["\$u_1/L\$" "\$-u_3/L\$"]
colors = [:blue,:orange]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$u_1/L, -u_3/L,\$", ylabel="\$F\$ [kip]", title="Displacements at point of force application")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo and Vu-Quoc (1986)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*F/(1e3), lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]]/L, [u1Ref[2,:],u3Ref[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/LeeFrameFollowerLoad_disp.pdf"))

println("Finished LeeFrameFollowerLoadPlotGenerator.jl")