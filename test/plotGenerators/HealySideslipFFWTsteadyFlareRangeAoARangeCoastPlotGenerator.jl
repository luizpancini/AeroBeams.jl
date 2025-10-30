using Plots, ColorSchemes

# Run the script
include("../examples/HealySideslipFFWTsteadyFlareRangeAoARangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTsteadyFlareRangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
ts = 12
fs = 16
lfs = 12
lw = 2
ms = 6
msw = 0
gr()

# Coast angle vs pitch angle for each flare angle
plt = plot(xlabel="Root pitch angle [deg]", ylabel="Coast angle [deg]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, xlims=[-20,30], ylims=[-90,90], yticks=-90:30:90, legend=:bottomright)
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Healy (2023) - Exp.")
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,Λ) in enumerate(ΛRange)
    plot!(θRange*180/π, -ϕHinge[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\Lambda = $(round(Int,Λ*180/π)) \\degree\$")
    if i==1
        plot!(flare10_num[1,:], flare10_num[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
        scatter!(flare10_exp[1,:], flare10_exp[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif i==2
        plot!(flare20_num[1,:], flare20_num[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
        scatter!(flare20_exp[1,:], flare20_exp[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    else
        plot!(flare30_num[1,:], flare30_num[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
        scatter!(flare30_exp[1,:], flare30_exp[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt)
savefig(string(absPath,"/HealySideslipFFWTsteadyFlareRangeAoARangeCoast_coastAngle.pdf"))

println("Finished HealySideslipFFWTsteadyFlareRangeAoARangeCoastPlotGenerator.jl")