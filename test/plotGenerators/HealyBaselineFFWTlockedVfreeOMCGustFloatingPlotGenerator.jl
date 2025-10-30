using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTlockedVfreeOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, 2, categorical=true)
ts = 12
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
gr()

# Root OOP bending moment increment
plt_ΔWRBM = plot(xlabel="\$\\omega_g\$ [Hz]", ylabel="ΔWRBM [N.m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, xlims=[0,15], legend=(0.25,0.65))
plot!([NaN],[NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label="Healy (2023) - Nastran")
plot!([NaN],[NaN], c=colors[1], lw=lw, shape=:circle, ms=ms, msw=msw, label="Free")
plot!([NaN],[NaN], c=colors[2], lw=lw, shape=:circle, ms=ms, msw=msw, label="Locked")
for (c,config) in enumerate(hingeConfigurations)
    ΔM2min = [ΔM2peaks[c,i][1] for i in 1:length(ωRange)]
    ΔM2max = [ΔM2peaks[c,i][2] for i in 1:length(ωRange)]
    plot!(ωRange, ΔM2min, c=colors[c], lw=lw, label=false)
    plot!(ωRange, ΔM2max, c=colors[c], lw=lw, label=false)
    if config == "locked"
        scatter!(ΔWRBM_locked_ref[1,:], ΔWRBM_locked_ref[2,:], c=colors[c], ms=ms, msw=msw, label=false)
    elseif config == "free"
        scatter!(ΔWRBM_free_ref[1,:], ΔWRBM_free_ref[2,:], c=colors[c], ms=ms, msw=msw, label=false)
    end
end
display(plt_ΔWRBM)
savefig(string(absPath,"/HealyBaselineFFWTlockedVfreeOMCGustFloating_DWRBM.pdf"))

println("Finished HealyBaselineFFWTlockedVfreeOMCGustFloatingPlotGenerator.jl")