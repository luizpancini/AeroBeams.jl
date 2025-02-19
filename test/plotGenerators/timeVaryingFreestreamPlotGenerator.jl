using Plots, ColorSchemes

# Run the script
include("../examples/timeVaryingFreestream.jl")

# Set paths
relPath = "/test/outputs/figures/timeVaryingFreestream"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λᵤRange)))
lw = 2
ms = 3
msw = 0
gr()

# Ratio of unsteady to quasi-steady cn over cycle
plt1 = plot(xlabel="\$t/T\$", ylabel="\$c_n/c_{n_{QS}}\$", xlims=[0,1], legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]/(2π*θ), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/timeVaryingFreestream_cn.pdf"))

# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt2)
savefig(string(absPath,"/timeVaryingFreestream_cm.pdf"))

# Relative wind acceleration over cycle
plt31 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt32 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/timeVaryingFreestream_Vdot.pdf"))

println("Finished timeVaryingFreestreamPlotGenerator.jl")