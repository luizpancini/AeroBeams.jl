using Plots, ColorSchemes

# Run the script
include("../examples/S10PazySteadyPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/S10PazySteadyPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
gr()

# Tip OOP displacement vs. airspeed
plt_tipOOP = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,70], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
scatter!([NaN], [NaN], ms=ms, msw=msw, c=:black, label="Exp.")
for (i,θ) in enumerate(θRange)
    plot!(URange, tipOOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label="\$\\theta = $(round(Int,θ*180/pi)) ^\\circ\$")
end
display(plt_tipOOP)
savefig(string(absPath,"/S10PazySteadyPitchRange_tipOOP.pdf"))

# Tip AoA vs. airspeed
plt_tipAOA = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[0,70], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipAoA[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/S10PazySteadyPitchRange_tipAoA.pdf"))

# Root axial strains vs. airspeed
plt_rootEps = plot(xlabel="Airspeed [m/s]", ylabel="Root axial strains (\$\\mu\$)", xlims=[0,70], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (i,θ) in enumerate(θRange)
    plot!(URange, rootEps[i,:]*1e6, c=colors[i], lw=lw, ls=:solid, label=false)
    if θ ≈ 0 + θoffset
        scatter!(steady_strains_aoa0_ref[1,:], steady_strains_aoa0_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 1*π/180 + θoffset
        scatter!(steady_strains_aoa1_ref[1,:], steady_strains_aoa1_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 3*π/180 + θoffset
        scatter!(steady_strains_aoa3_ref[1,:], steady_strains_aoa3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 5*π/180 + θoffset
        scatter!(steady_strains_aoa5_ref[1,:], steady_strains_aoa5_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 7*π/180 + θoffset
        scatter!(steady_strains_aoa7_ref[1,:], steady_strains_aoa7_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 10*π/180 + θoffset
        scatter!(steady_strains_aoa10_ref[1,:], steady_strains_aoa10_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_rootEps)
savefig(string(absPath,"/S10PazySteadyPitchRange_rootEps.pdf"))

println("Finished S10PazySteadyPitchRangePlotGenerator.jl")