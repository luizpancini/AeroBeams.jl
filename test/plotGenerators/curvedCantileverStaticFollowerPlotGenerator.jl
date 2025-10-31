using Plots

# Run the script
include("../examples/curvedCantileverStaticFollower.jl")

# Set paths
relPath = "/test/outputs/figures/curvedCantileverStaticFollower"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,σ2plot=vcat(0.2:0.2:1),legendPos=(0.5,0.6),view=(135,15),save=true,savePath=string(relPath,"/curvedCantileverStaticFollower_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = [:blue,:green,:orange]
labels = ["\$-u_1\$" "\$u_2\$" "\$u_3\$"]
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
gr()

# Plot normalized tip displacements over load steps
plt1 = plot(xlabel="\$F\$ [lb]", ylabel="Tip displacements [in]", xlims=[0,3000], xticks=collect(0:500:F), yticks=collect(-60:20:80), legend_position=(0.45, 0.4), tickfont=font(ts), guidefont=font(fs), legendfontsize=10)
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1986)")
for i=1:3
    plot!([NaN], [NaN], lc=colors[i], m=colors[i],  lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (y, c) in zip([-tip_u1, tip_u2, tip_u3], colors)
    plot!(σVector*F, y, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_ref[1,:],u2_ref[1,:],u3_ref[1,:]], [u1_ref[2,:],u2_ref[2,:],u3_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/curvedCantileverStaticFollower_disp.pdf"))

println("Finished curvedCantileverStaticFollower.jl")