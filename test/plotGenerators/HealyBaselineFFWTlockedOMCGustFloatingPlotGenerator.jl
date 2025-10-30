using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTlockedOMCGustFloating.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTlockedOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(ωRange), categorical=true)
ts = 10
fs = 16
lw = 2
gr()

# Root OOP bending moment increment
plt_ΔWRBM = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
for (i,ω) in enumerate(ωRange)
    plot!(t[i], -(M2root[i] .- M2root[i][1]), c=colors[i], lw=lw, ls=:solid, label="\$\\omega_g = $ω\$ Hz")
    if ω==1
        plot!(ΔWRBM_Healy_1Hz[1,:],ΔWRBM_Healy_1Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif ω==3
        plot!(ΔWRBM_Healy_3Hz[1,:],ΔWRBM_Healy_3Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif ω==10
        plot!(ΔWRBM_Healy_10Hz[1,:],ΔWRBM_Healy_10Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    end
end
display(plt_ΔWRBM)
savefig(string(absPath,"/HealyBaselineFFWTlockedOMCGustFloating_DWRBM.pdf"))

println("Finished HealyBaselineFFWTlockedOMCGustFloatingPlotGenerator.jl")