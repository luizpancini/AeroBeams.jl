using Plots, ColorSchemes

# Run the script
include("../examples/HealyLCOFFWTfreeFlutterAoARange.jl")

# Set paths
relPath = "/test/outputs/figures/HealyLCOFFWTfreeFlutterAoARange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
ms = 3
msw = 0
gr()

# Fold angle
plt_fold = plot(xlabel="Airspeed [m/s]", ylabel="Fold angle [deg]", xlims=[0,30], ylims=[-90,30], yticks=-90:30:30)
# plot!([NaN],[NaN], lc=:black, ls=:dashdot, lw=lw, label="Healy (2023)")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (j,θ) in enumerate(θRange)
    plot!(URange, ϕHinge[j,:], c=colors[j], lw=lw, label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
end
display(plt_fold)
savefig(string(absPath,"/HealyLCOFFWTfreeFlutterAoARange_fold.pdf"))

# V-g-f - first 2 modes
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,30], ylims=[0,3], legend=:bottomright)
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,30], ylims=[-0.5,0.25], legend=:bottomleft)
# plot!(plt_Vf, [NaN],[NaN], c=:black, ls=:dashdot, lw=lw, label="Healy (2023)")
plot!(plt_Vf, [NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (j,θ) in enumerate(θRange)
    for mode in 1:2
        plot!(plt_Vf, URange, modeFrequencies[j,mode]/(2*π), c=colors[j], lw=lw, label=false)
        plot!(plt_Vg, URange, modeDampingRatios[j,mode], c=colors[j], lw=lw, label=false)
    end
    plot!(plt_Vg,[NaN],[NaN], c=colors[j], ls=:solid, lw=lw, label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/HealyLCOFFWTfreeFlutterAoARange_Vgf_first2modes.pdf"))

# # V-g-f - first 5 modes
# plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,30], ylims=[0,40], legend=:bottomright)
# plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,30], ylims=[-0.5,0.25], legend=:bottomleft)
# # plot!(plt_Vf, [NaN],[NaN], c=:black, ls=:dashdot, lw=lw, label="Healy (2023)")
# plot!(plt_Vf, [NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
# for (j,θ) in enumerate(θRange)
#     for mode in 1:8
#         plot!(plt_Vf, URange, modeFrequencies[j,mode]/(2*π), c=colors[j], lw=lw, label=false)
#         plot!(plt_Vg, URange, modeDampingRatios[j,mode], c=colors[j], lw=lw, label=false)
#     end
#     plot!(plt_Vg,[NaN],[NaN], c=colors[j], ls=:solid, lw=lw, label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
# end
# plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
# display(plt_Vgf)
# savefig(string(absPath,"/HealyLCOFFWTfreeFlutterAoARange_Vgf_first5modes.pdf"))

println("Finished HealyLCOFFWTfreeFlutterAoARangePlotGenerator.jl")