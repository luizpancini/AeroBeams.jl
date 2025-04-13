using Plots, ColorSchemes

# Run the script
include("../examples/HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
plotTitle = false
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lw = 2
ms = 6
msw = 0
gr()

# Coast angle vs sideslip angle for each AoA, for several wingtip twist angles
for (i,φ) in enumerate(φRange)
    if i==1
        figurePath = "/HealyFFWTsteadySideslipRangeAoARangeCoast_phim6.pdf"
        aoam3_exp = phim6_aoam3_exp
        aoam3_num = phim6_aoam3_num
        aoa0_exp = phim6_aoa0_exp
        aoa0_num = phim6_aoa0_num
        aoa3_exp = phim6_aoa3_exp
        aoa3_num = phim6_aoa3_num
        aoa6_exp = phim6_aoa6_exp
        aoa6_num = phim6_aoa6_num
        aoa9_exp = phim6_aoa9_exp
        aoa9_num = phim6_aoa9_num
        aoa12_exp = phim6_aoa12_exp
        aoa12_num = phim6_aoa12_num
    elseif i==2
        figurePath = "/HealyFFWTsteadySideslipRangeAoARangeCoast_phi0.pdf"
        aoam3_exp = phi0_aoam3_exp
        aoam3_num = phi0_aoam3_num
        aoa0_exp = phi0_aoa0_exp
        aoa0_num = phi0_aoa0_num
        aoa3_exp = phi0_aoa3_exp
        aoa3_num = phi0_aoa3_num
        aoa6_exp = phi0_aoa6_exp
        aoa6_num = phi0_aoa6_num
        aoa9_exp = phi0_aoa9_exp
        aoa9_num = phi0_aoa9_num
        aoa12_exp = phi0_aoa12_exp
        aoa12_num = phi0_aoa12_num
    else
        figurePath = "/HealyFFWTsteadySideslipRangeAoARangeCoast_phi9.pdf"
        aoam3_exp = phi9_aoam3_exp
        aoam3_num = phi9_aoam3_num
        aoa0_exp = phi9_aoa0_exp
        aoa0_num = phi9_aoa0_num
        aoa3_exp = phi9_aoa3_exp
        aoa3_num = phi9_aoa3_num
        aoa6_exp = phi9_aoa6_exp
        aoa6_num = phi9_aoa6_num
        aoa9_exp = phi9_aoa9_exp
        aoa9_num = phi9_aoa9_num
        aoa12_exp = phi9_aoa12_exp
        aoa12_num = phi9_aoa12_num
    end
    plt = plot(xlabel="Sideslip angle [deg]", ylabel="Coast angle [deg]", xlims=[-10,30], ylims=[-150,150], yticks=-150:30:150, tickfont=font(ts), guidefont=font(fs))
    if i==1
        plot!(legendfontsize=7, legend=:topleft)
    else
        plot!(legend=false)
    end
    if plotTitle
        plot!(title="Wingtip twist = \$ $(round(Int,φ*180/π)) \\degree\$")
    end
    scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Exp. - Healy (2023)")
    plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Num. - Healy (2023)")
    plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
    for (j,θ) in enumerate(θRange)
        plot!([NaN], [NaN], lw=lw, ls=:solid, ms=ms, msw=msw, c=colors[j], label="\$\\theta = $(round(Int,θ*180/π)) \\degree\$")
    end
    for (j,θ) in enumerate(θRange)
        plot!(βRange*180/π, -ϕHinge[i,j,:], lw=lw, ls=:solid, c=colors[j], label=false)
        if j==1
            plot!(aoam3_num[1,:], aoam3_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoam3_exp[1,:], aoam3_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        elseif j==2
            plot!(aoa0_num[1,:], aoa0_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa0_exp[1,:], aoa0_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        elseif j==3
            plot!(aoa3_num[1,:], aoa3_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa3_exp[1,:], aoa3_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        elseif j==4
            plot!(aoa6_num[1,:], aoa6_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa6_exp[1,:], aoa6_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        elseif j==5
            plot!(aoa9_num[1,:], aoa9_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa9_exp[1,:], aoa9_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        else
            plot!(aoa12_num[1,:], aoa12_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa12_exp[1,:], aoa12_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)     
        end
    end
    display(plt)
    savefig(string(absPath,figurePath))
end

println("Finished HealyFFWTsteadySideslipRangeAoARangeCoastPlotGenerator.jl")