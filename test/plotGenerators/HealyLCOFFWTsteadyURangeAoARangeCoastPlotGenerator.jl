using Plots, ColorSchemes

# Run the script
include("../examples/HealyLCOFFWTsteadyURangeAoARangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/HealyLCOFFWTsteadyURangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 12
fs = 16
lfs = 9
lw = 2
ms = 6
msw = 0
gr()
lstyle = [:dash, :solid]

# Select airspeed indices for plots
kPlot = [1, length(URange)]

# Plots
plt_OOP = plot(xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized OOP position", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomleft)
plot!([NaN], [NaN], lw=lw, ls=lstyle[1], c=:black, label="\$U = $(URange[kPlot[1]])\$ m/s")
plot!([NaN], [NaN], lw=lw, ls=lstyle[2], c=:black, label="\$U = $(URange[kPlot[2]])\$ m/s")
for (j,θ) in enumerate(θRange)
    plot!([NaN], [NaN], lw=lw, c=colors[j], label="\$\\theta = $(round(θ*180/π,digits=1)) ^\\circ\$")
end
plt_p1 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_1\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_p2 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_2\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_p3 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_3\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_M2 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_AOA = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_cn = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)

for (j,θ) in enumerate(θRange)
    for (kk, k) in enumerate(kPlot)
        U = URange[k]
        # OOP displacement
        plot!(plt_OOP, (r_n1.+u1[j,k])/x1_n[end], (r_n3.+u3[j,k])/x1_n[end], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # p1
        plot!(plt_p1, x1_n/x1_n[end], p1[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # p2
        plot!(plt_p2, x1_n/x1_n[end], p2[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # p3
        plot!(plt_p3, x1_n/x1_n[end], p3[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # Bending moment
        plot!(plt_M2, x1_n/x1_n[end], M2[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # Angle of attack
        plot!(plt_AOA, x1_e/x1_n[end], α[j,k]*180/π, lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # cn
        plot!(plt_cn, x1_e/x1_n[end], cn[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
    end
end

# Display and save figures
display(plt_OOP)
display(plt_p1)
display(plt_p2)
display(plt_p3)
display(plt_M2)
display(plt_AOA)
display(plt_cn)
savefig(plt_OOP, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_u3.pdf"))
savefig(plt_p1, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_p1.pdf"))
savefig(plt_p2, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_p2.pdf"))
savefig(plt_p3, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_p3.pdf"))
savefig(plt_M2, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_M2.pdf"))
savefig(plt_AOA, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_aoa.pdf"))
savefig(plt_cn, string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_cn.pdf"))

# Coast angle vs airspeed for each AoA
plt = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", xlims=[0,30], ylims=[-90,60], yticks=-90:15:60, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topleft)
plot!([NaN],[NaN], c=:black, ls=:solid, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label="Healy (2023) - Exp.")
plot!([NaN],[NaN], c=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
for (j,θ) in enumerate(θRange)
    plot!(URange, -ϕHinge[j,:], lw=lw, ls=:solid, c=colors[j], label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
    if θ == 2.5*π/180
        scatter!(aoa_25_exp[1,:], aoa_25_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_25_num[1,:], aoa_25_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    elseif θ == 5.0*π/180
        scatter!(aoa_50_exp[1,:], aoa_50_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_50_num[1,:], aoa_50_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    elseif θ == 7.5*π/180
        scatter!(aoa_75_exp[1,:], aoa_75_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_75_num[1,:], aoa_75_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    elseif θ == 10*π/180
        scatter!(aoa_10_exp[1,:], aoa_10_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_10_num[1,:], aoa_10_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    end
end
display(plt)
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_coastAngle.pdf"))

println("Finished HealyLCOFFWTsteadyURangeAoARangeCoastPlotGenerator.jl")