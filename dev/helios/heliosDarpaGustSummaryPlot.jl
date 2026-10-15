
# Gust intensity range
γRange = [0.1, 0.15, 0.20, 0.25]

# Instability gust duration values
τi_AF = [27.5, 9.0, 6.5, 5.0]
τi_DS = [10.0, 5.0, 3.5, 3.0]

# Airspeed [ft/s] and payload [lb]
U = 40
P = 150

# Set paths
relPath = "/dev/helios/figures/heliosDarpaGustSummaryPlot"
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

plt_Hinst = plot(xlabel="\$\\sigma_g\$ [%]", ylabel="\$H\$ for instability [semichords]", xlims=[9.95,25.05], ylims=[0,300], xticks=vcat(10:5:25), tickfont=font(ts), guidefont=font(fs), legend=:topright, legendfontsize=lfs)
plot!(γRange*100, τi_AF*U/4, marker=:circle, c=colors[1], lw=lw, ms=ms, msw=msw, label="Attached flow")
plot!(γRange*100, τi_DS*U/4, marker=:circle, c=colors[2], lw=lw, ms=ms, msw=msw, label="Dynamic stall")
display(plt_Hinst)
savefig(plt_Hinst,string(absPath,"/heliosDarpaGustSummaryPlot_P",round(P),"_U",round(U),".pdf"))
