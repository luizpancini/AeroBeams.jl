using Plots

# Run the script
include("../examples/timeVaryingFreestreamAndPitch.jl")

# Set paths
relPath = "/test/outputs/figures/timeVaryingFreestreamAndPitch"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,showScale=false,plotAeroSurf=false,plotLimits=[(0,L),(-L/2,L/2),(-L/2,L/2)],save=true,savePath=string(relPath,"/timeVaryingFreestreamAndPitch_deformation.gif"),displayProgress=true)

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
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]./cn_qs(t[i][rangeLastCycle[i]]), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_cn.pdf"))

# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt2)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_cm.pdf"))

# Relative wind velocity over cycle
plt31 = plot(xlabel="\$t/T\$", ylabel="\$V_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], V2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], V2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt32 = plot(xlabel="\$t/T\$", ylabel="\$V_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], V3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], V3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt33 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Ω1[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ω1Analytical.(t[i][rangeLastCycle[i][1:10:end]]), c=colors[i], ms=ms, msw=msw, label=false)
end
plt3 = plot(plt31,plt32,plt33, layout=(3,1))
display(plt3)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_velocities.pdf"))

# Relative wind acceleration over cycle
plt41 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt42 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt43 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1Analytical.(t[i][rangeLastCycle[i][1:10:end]]), c=colors[i], ms=ms, msw=msw, label=false)
end
plt4 = plot(plt41,plt42,plt43, layout=(3,1))
display(plt4)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_accelerations.pdf"))

# Effective angle of attack over cycle
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\alpha_E\$ [deg]", xlims=[0,1], legend=:bottomleft)
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], αₑ[i][rangeLastCycle[i]]*180/π, c=colors[i], lw=lw, label="λ = $λᵤ")
end
display(plt5)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_alpha.pdf"))

println("Finished timeVaryingFreestreamAndPitchPlotGenerator.jl")