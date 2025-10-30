using Plots, ColorSchemes

# Run the script
include("../examples/rotorDStest.jl")

# Set paths
relPath = "/test/outputs/figures/rotorDStest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(dynProblem,refBasis="I",plotFrequency=round(Int,τ/40/Δt),fps=30,view=(30,30),plotLimits=([-L,L],[-L,L],[-L,L]),plotDistLoads=true,save=true,savePath=string(relPath,"/rotorDStest_deformation.gif"),displayProgress=true)

# Plot configurations
colors = cgrad(:rainbow, length(elemRangePlot), categorical=true)
labels = ["Root element", "Element $(elemRangePlot[2])", "Tip element"]
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
gr()

# Indices defining time range to plot
lastCycles2plot = nCycles
ind = findfirst(x-> x>=t[end]-lastCycles2plot*τ, t)
range2plot = ind:length(t)

# Prescribed root pitch angle
θprescribed = (a₀-a₁.+θ.(t))*180/π
θprescribedRange2plot = 1:50:length(t)

# Pitch angle
plt_alpha = plot(xlabel="Cycle", ylabel="Pitch angle [deg]", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
scatter!(t[θprescribedRange2plot]/τ, θprescribed[θprescribedRange2plot], c=:black, ms=ms, label="Prescribed at the root")
for (i,e) in enumerate(elemRangePlot)
    plot!(t[range2plot]/τ, α[i][range2plot]*180/π, c=colors[i], lw=lw, label=labels[i])
end
display(plt_alpha)
savefig(string(absPath,"/rotorDStest_alpha.pdf"))

# Relative airspeed
plt_U = plot(xlabel="Cycle", ylabel="Airspeed [m/s]", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t[range2plot]/τ, U*ones(length(range2plot)), c=:black, ls=:dash, lw=lw, label="Reference airspeed")
for (i,e) in enumerate(elemRangePlot)
    plot!(t[range2plot]/τ, U∞rel[i][range2plot], c=colors[i], lw=lw, label=labels[i])
end
display(plt_U)
savefig(string(absPath,"/rotorDStest_U.pdf"))

# Vertical displacement
plt_u3 = plot(xlabel="Cycle", ylabel="\$u_3\$ [m]", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t[range2plot]/τ, u3[i][range2plot], c=colors[i], lw=lw, label=false)
end
display(plt_u3)
savefig(string(absPath,"/rotorDStest_u3.pdf"))

# cn vs time
plt_cnt = plot(xlabel="Cycle", ylabel="\$c_n\$", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t[range2plot]/τ, cn[i][range2plot], c=colors[i], lw=lw, label=false)
end
display(plt_cnt)
savefig(string(absPath,"/rotorDStest_cnt.pdf"))

# cm vs time
plt_cmt = plot(xlabel="Cycle", ylabel="\$c_m\$", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t[range2plot]/τ, cm[i][range2plot], c=colors[i], lw=lw, label=false)
end
display(plt_cmt)
savefig(string(absPath,"/rotorDStest_cmt.pdf"))

# ct vs time
plt_ctt = plot(xlabel="Cycle", ylabel="\$c_t\$", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,e) in enumerate(elemRangePlot)
    plot!(t[range2plot]/τ, ct[i][range2plot], c=colors[i], lw=lw, label=false)
end
display(plt_ctt)
savefig(string(absPath,"/rotorDStest_ctt.pdf"))

# cl vs α
plt_cla = plot(xlabel="Pitch angle [deg]", ylabel="\$c_l\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
scatter!(clRef[1,:], clRef[2,:], c=:black, ms=ms, label="Exp. - McAlister et al (1982)")
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i][range2plot]*180/π, cl[i][range2plot], c=colors[i], lw=lw, label=labels[i])
end
display(plt_cla)
savefig(string(absPath,"/rotorDStest_cla.pdf"))

# cm vs α
plt_cma = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legend=:bottomleft)
scatter!(cmRef[1,:], cmRef[2,:], c=:black, ms=ms, label=false)
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i][range2plot]*180/π, cm[i][range2plot], c=colors[i], lw=lw, label=false)
end
display(plt_cma)
savefig(string(absPath,"/rotorDStest_cma.pdf"))

# cd vs α
plt_cda = plot(xlabel="Pitch angle [deg]", ylabel="\$c_d\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=12)
scatter!(cdRef[1,:], cdRef[2,:], c=:black, ms=ms, label=false)
for (i,e) in enumerate(elemRangePlot)
    plot!(α[i][range2plot]*180/π, cdrag[i][range2plot], c=colors[i], lw=lw, label=false)
end
display(plt_cda)
savefig(string(absPath,"/rotorDStest_cda.pdf"))

println("Finished rotorDStestPlotGenerator.jl")