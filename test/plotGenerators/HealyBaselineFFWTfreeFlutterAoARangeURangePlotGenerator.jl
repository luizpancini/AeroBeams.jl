using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTfreeFlutterAoARangeURange.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTfreeFlutterAoARangeURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
ms = 3
msw = 0
gr()

# Fold angle - without tip loss
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Fold angle [deg]", xlims=[0,40], ylims=[-90,30], yticks=-90:30:30, title="Without tip loss")
plot!([NaN],[NaN], lc=:black, ls=:dashdot, lw=lw, label="Healy (2023) - ST")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (j,θ) in enumerate(θRange)
    plot!(fold_ST[2*j-1,:],fold_ST[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(URange, ϕHinge[2,j,:], c=colors[j], lw=lw, label="\$\\theta = $(round(Int,θ*180/π)) \\degree\$")
end
display(plt1)
savefig(string(absPath,"/HealyBaselineFFWTfreeFlutterAoARangeURange_fold_tipLoss0.pdf"))

# Fold angle - with tip loss
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Fold angle [deg]", xlims=[0,40], ylims=[-90,30], yticks=-90:30:30, title="With tip loss")
plot!([NaN],[NaN], lc=:black, ls=:dashdot, lw=lw, label="Healy (2023) - VLM")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (j,θ) in enumerate(θRange)
    plot!(fold_VLM[2*j-1,:],fold_VLM[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(URange, ϕHinge[1,j,:], c=colors[j], lw=lw, label="\$\\theta = $(round(Int,θ*180/π)) \\degree\$")
end
display(plt2)
savefig(string(absPath,"/HealyBaselineFFWTfreeFlutterAoARangeURange_fold_tipLoss1.pdf"))

# V-g-f - without tip loss
range2plot = 1:findfirst(x -> x >= 22, URange)
plt31 = plot(ylabel="Frequency [Hz]", xlims=[0,30], ylims=[0,5], title="Without tip loss", legend=:bottomright)
plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,30], ylims=[-0.5,0.25], legend=:bottomleft)
plot!(plt31, [NaN],[NaN], lc=:black, ls=:dashdot, lw=lw, label="Healy (2023) - ST")
plot!(plt31, [NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (j,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(plt31, URange, modeFrequencies[2,j,mode]/(2*π), c=colors[j], lw=lw, label=false)
        plot!(plt32, URange, modeDampingRatios[2,j,mode], c=colors[j], lw=lw, label=false)
    end
    plot!(plt31, freq1_ST[2*j-1,:], freq1_ST[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt31, freq2_ST[2*j-1,:], freq2_ST[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt32, damp1_ST[2*j-1,:], damp1_ST[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt32, damp2_ST[2*j-1,:], damp2_ST[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt32,[NaN],[NaN], c=colors[j], ls=:solid, lw=lw, label="\$\\theta = $(round(Int,θ*180/π)) \\degree\$")
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/HealyBaselineFFWTfreeFlutterAoARangeURange_Vgf_tipLoss0.pdf"))


# V-g-f - with tip loss
plt41 = plot(ylabel="Frequency [Hz]", xlims=[0,30], ylims=[0,5], title="With tip loss", legend=:bottomright)
plt42 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,30], ylims=[-0.5,0.25], legend=:bottomleft)
plot!(plt41, [NaN],[NaN], lc=:black, ls=:dashdot, lw=lw, label="Healy (2023) - VLM")
plot!(plt41, [NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (j,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(plt41, URange, modeFrequencies[1,j,mode]/(2*π), c=colors[j], lw=lw, label=false)
        plot!(plt42, URange, modeDampingRatios[1,j,mode], c=colors[j], lw=lw, label=false)
    end
    plot!(plt41, freq1_VLM[2*j-1,:], freq1_VLM[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt41, freq2_VLM[2*j-1,:], freq2_VLM[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt42, damp1_VLM[2*j-1,:], damp1_VLM[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt42, damp2_VLM[2*j-1,:], damp2_VLM[2*j,:], c=colors[j], ls=:dashdot, lw=lw, label=false)
    plot!(plt42,[NaN],[NaN], c=colors[j], ls=:solid, lw=lw, label="\$\\theta = $(round(Int,θ*180/π)) \\degree\$")
end
plt4 = plot(plt41,plt42, layout=(2,1))
display(plt4)
savefig(string(absPath,"/HealyBaselineFFWTfreeFlutterAoARangeURange_Vgf_tipLoss1.pdf"))

println("Finished HealyBaselineFFWTfreeFlutterAoARangeURangePlotGenerator.jl")