using Plots, ColorSchemes

# Run the script
include("../examples/HealyFFWTsteadyFlareRangeURangeAoARangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/HealyFFWTsteadyFlareRangeURangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
ms = 6
msw = 0
gr()

# Coast angle vs airspeed for each AoA, for several flare angles
for (i,Λ) in enumerate(ΛRange)
    if i==1
        figurePath = "/HealyFFWTsteadyURangeAoARangeCoast_phi_flare10.pdf"
        aoa0_exp = flare10_aoa0_exp
        aoa0_num = flare10_aoa0_num
        aoa5_exp = flare10_aoa5_exp
        aoa5_num = flare10_aoa5_num
        aoa10_exp = flare10_aoa10_exp
        aoa10_num = flare10_aoa10_num
    elseif i==2
        figurePath = "/HealyFFWTsteadyURangeAoARangeCoast_phi_flare20.pdf"
        aoa0_exp = flare20_aoa0_exp
        aoa0_num = flare20_aoa0_num
        aoa5_exp = flare20_aoa5_exp
        aoa5_num = flare20_aoa5_num
        aoa10_exp = flare20_aoa10_exp
        aoa10_num = flare20_aoa10_num
    else
        figurePath = "/HealyFFWTsteadyURangeAoARangeCoast_phi_flare30.pdf"
        aoa0_exp = flare30_aoa0_exp
        aoa0_num = flare30_aoa0_num
        aoa5_exp = flare30_aoa5_exp
        aoa5_num = flare30_aoa5_num
        aoa10_exp = flare30_aoa10_exp
        aoa10_num = flare30_aoa10_num
    end
    plt = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", xlims=[10,35], ylims=[-90,45], yticks=-90:30:45, legendfontsize=10, legend=:bottomright)
    scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Healy (2023) - Exp.")
    plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
    plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
    for (j,θ) in enumerate(θRange)
        plot!(URange, -ϕHinge[i,j,:], lw=lw, ls=:solid, c=colors[j], label="\$\\theta = $(round(Int,θ*180/π)) \\degree\$")
        if j==1
            plot!(aoa0_num[1,:], aoa0_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa0_exp[1,:], aoa0_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        elseif j==2
            plot!(aoa5_num[1,:], aoa5_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa5_exp[1,:], aoa5_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        else
            plot!(aoa10_num[1,:], aoa10_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa10_exp[1,:], aoa10_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        end
    end
    display(plt)
    savefig(string(absPath,figurePath))
end

println("Finished HealyFFWTsteadyURangeAoARangeCoastPlotGenerator.jl")