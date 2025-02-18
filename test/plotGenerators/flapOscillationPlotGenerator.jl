using Plots

# Run the script
include("../examples/flapOscillation.jl")

# Set paths
relPath = "/test/outputs/figures/flapOscillation"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
rangeLastCycle = findfirst(x->x==tf-T,t):length(t)
ts = 10
fs = 16
lw = 2
ms = 3
gr()

# cn vs δ
plt1 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$c_n/\\pi\$", xlims=[-3,3], ylims=[-0.075,0.075], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!(δ.(t[rangeLastCycle])*180/π, cn[rangeLastCycle]/π, c=:black, lw=lw, label="AeroBeams")
plot!(cnRefMod[1,:], cnRefMod[2,:]/π, c=:black, ls=:dashdot, lw=lw, label="Model - Leishman (2006)")
scatter!(cnExp[1,:], cnExp[2,:]/π, c=:black, ms=ms, msw=0, label="Exp. - Tijdeman & Schippers (1973)")
display(plt1)
savefig(string(absPath,"/flapOscillation_cn.pdf"))

# cm vs δ
plt2 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$-2c_m/\\pi\$", xlims=[-3,3], ylims=[-0.03,0.03], tickfont=font(ts), guidefont=font(fs))
plot!(δ.(t[rangeLastCycle])*180/π, -2*cm[rangeLastCycle]/π, c=:black, lw=lw, label=false)
plot!(cmRefMod[1,:], -2*cmRefMod[2,:]/π, c=:black, ls=:dashdot, lw=lw, label=false)
scatter!(cmExp[1,:], -2*cmExp[2,:]/π, c=:black, ms=ms, msw=0, label=false)
display(plt2)
savefig(string(absPath,"/flapOscillation_cm.pdf"))

println("Finished flapOscillation.jl")