using Plots, Measures

# Payload range
PRange = vcat(0:25:225)

# Instability flap deflection increments
δi_AF_p =  [28, 27, 27, 26, 25, 25, 24, 23, 22, 0]
δi_DS_p =  [13, 12, 12, 12, 11, 11, 10,  9,  3, 0]
δi_AF_n = -[17, 17, 18, 18, 19, 19, 20, 21, 22, 0]
δi_DS_n = -[ 9,  9,  9,  9,  9,  9,  9,  8,  3, 0]

# Airspeed [ft/s]
U = 40

# Set paths
relPath = "/dev/helios/figures/heliosPitchManeuverSummaryPlot"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 5
msw = 0
colors = cgrad(:rainbow, 2, categorical=true)

plt_Hinst = plot(xlabel="Payload [lb]", ylabel="\$\\Delta\\delta_f\$ for instability [deg]", xlims=[0,225], ylims=[-30,30], xticks=vcat(0:25:225), yticks=vcat(-30:5:30), tickfont=font(ts), guidefont=font(fs), legend=(0.1,0.12), legendfontsize=lfs, left_margin=2mm, right_margin=2mm)
plot!(PRange, δi_AF_p, marker=:circle, c=colors[1], lw=lw, ms=ms, msw=msw, label="Attached flow")
plot!(PRange, δi_DS_p, marker=:circle, c=colors[2], lw=lw, ms=ms, msw=msw, label="Dynamic stall")
plot!(PRange, δi_AF_n, marker=:circle, c=colors[1], lw=lw, ms=ms, msw=msw, label=false)
plot!(PRange, δi_DS_n, marker=:circle, c=colors[2], lw=lw, ms=ms, msw=msw, label=false)
display(plt_Hinst)
savefig(plt_Hinst,string(absPath,"/heliosPitchManeuverSummaryPlot_U",round(U),".pdf"))
