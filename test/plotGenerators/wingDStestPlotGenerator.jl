using Plots, ColorSchemes

# Run the script
include("../examples/wingDStest.jl")

# Set paths
relPath = "/test/outputs/figures/wingDStest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animations
for (i,aeroSolver) in enumerate(aeroSolvers)
    anim = plot_dynamic_deformation(dynProblem[i],refBasis="I",followAssembly=true,scale=1,plotFrequency=round(Int,τ/50/Δt),view=(60,15),plotLimits=([-L/4,L],[-L/4,L],[-L/4,L]),plotDistLoads=true,save=true,savePath=string(relPath,"/wingDStest_deformation_",aeroSolver.name,".gif"),displayProgress=true)
    display(anim)
end

# Plot configurations
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)
elemLabels = ["Root", "Tip"]
lstyles = [:solid, :dash]
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
fps = 30
DPI = 300
gr()

# Number of cycles to plot
lastCycles2plot = 1

# Pitch angle
plt_alpha = plot(xlabel="Cycle", ylabel="Pitch angle [deg]", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
scatter!([NaN], [NaN], c=:black, ms=ms, label="Prescribed")
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Prescribed root pitch angle
    θprescribed = (a₀-a₁.+θ.(t[i]))*180/π
    θprescribedRange2plot = 1:50:length(t[i])
    scatter!(t[i][θprescribedRange2plot]/τ, θprescribed[θprescribedRange2plot], c=:black, ms=ms, label=false)
    #
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(t[i][range2plot]/τ .- (nCycles-lastCycles2plot), α[i][j][range2plot]*180/π, c=colors[i], lw=lw, ls=lstyles[j], label=string(aeroSolver.name," - ",elemLabels[j]))
    end
end
display(plt_alpha)
savefig(string(absPath,"/wingDStest_alpha.pdf"))

# cn vs time
plt_cnt = plot(xlabel="Cycle", ylabel="\$c_n\$", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(t[i][range2plot]/τ .- (nCycles-lastCycles2plot), cn[i][j][range2plot], c=colors[i], lw=lw, ls=lstyles[j], label=string(aeroSolver.name," - ",elemLabels[j]))
    end
end
display(plt_cnt)
savefig(string(absPath,"/wingDStest_cnt.pdf"))

# cm vs time
plt_cmt = plot(xlabel="Cycle", ylabel="\$c_m\$", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(t[i][range2plot]/τ .- (nCycles-lastCycles2plot), cm[i][j][range2plot], c=colors[i], lw=lw, ls=lstyles[j], label=false)
    end
end
display(plt_cmt)
savefig(string(absPath,"/wingDStest_cmt.pdf"))

# ct vs time
plt_ctt = plot(xlabel="Cycle", ylabel="\$c_t\$", xlims=[0,lastCycles2plot], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(t[i][range2plot]/τ .- (nCycles-lastCycles2plot), ct[i][j][range2plot], c=colors[i], lw=lw, ls=lstyles[j], label=false)
    end
end
display(plt_ctt)
savefig(string(absPath,"/wingDStest_ctt.pdf"))

# cl vs α for root element only
labels = ["Attached-flow" " " "Dynamic stall"]
plt_cla = plot(xlabel="Pitch angle [deg]", ylabel="\$c_l\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs, dpi=DPI)
scatter!(clRef[1,:], clRef[2,:], c=:black, ms=ms, label="Exp. - McAlister et al (1982)")
for (i,aeroSolver) in enumerate(aeroSolvers)
    if i==2
        continue
    end
    plot!([NaN], [NaN], c=colors[i], lw=lw, label=labels[i])
end
plt_cla_base = deepcopy(plt_cla)
for (i,aeroSolver) in enumerate(aeroSolvers)
    if i==2
        continue
    end
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    plot!(α[i][1][range2plot]*180/π, cl[i][1][range2plot], c=colors[i], lw=lw, ls=lstyles[1], label=false)
end
display(plt_cla)
savefig(string(absPath,"/wingDStest_cla_root.pdf"))

# cl vs alpha animations or root element only
plt_cl_anim = plt_cla_base
ind = findfirst(x-> x>=t[1][end]-lastCycles2plot*τ, t[1])
N = length(t[1])
tCycle = t[1][ind:N]
plotFrequency = 10
anim = @animate for (k,timeNow) in enumerate(tCycle)
    if k > 1 && rem(k,plotFrequency) > 0
        continue
    end
    plot(plt_cl_anim)
    for (i,aeroSolver) in enumerate(aeroSolvers)
        if i==2
            continue
        end
        plot!(α[i][1][ind:ind+k-1]*180/π, cl[i][1][ind:ind+k-1], c=colors[i], lw=lw, label=false)
    end
end
gif_handle = gif(anim, string(absPath,"/wingDStest_cla_root.gif"), fps=fps)
display(gif_handle)

# cm vs α for root element only
plt_cma = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs, dpi=DPI)
scatter!(cmRef[1,:], cmRef[2,:], c=:black, ms=ms, label=false)
plt_cma_base = deepcopy(plt_cma)
for (i,aeroSolver) in enumerate(aeroSolvers)
    if i==2
        continue
    end
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    plot!(α[i][1][range2plot]*180/π, cm[i][1][range2plot], c=colors[i], lw=lw, ls=lstyles[1], label=false)
end
display(plt_cma)
savefig(string(absPath,"/wingDStest_cma_root.pdf"))

# cm vs alpha animations or root element only
plt_cm_anim = plt_cma_base
ind = findfirst(x-> x>=t[1][end]-lastCycles2plot*τ, t[1])
N = length(t[1])
tCycle = t[1][ind:N]
plotFrequency = 10
anim = @animate for (k,timeNow) in enumerate(tCycle)
    if k > 1 && rem(k,plotFrequency) > 0
        continue
    end
    plot(plt_cm_anim)
    for (i,aeroSolver) in enumerate(aeroSolvers)
        if i==2
            continue
        end
        plot!(α[i][1][ind:ind+k-1]*180/π, cm[i][1][ind:ind+k-1], c=colors[i], lw=lw, label=false)
    end
end
gif_handle = gif(anim, string(absPath,"/wingDStest_cma_root.gif"), fps=fps)
display(gif_handle)

# cl vs α
plt_cla = plot(xlabel="Pitch angle [deg]", ylabel="\$c_l\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
scatter!(clRef[1,:], clRef[2,:], c=:black, ms=ms, label="Exp. - McAlister et al (1982)")
for (i,aeroSolver) in enumerate(aeroSolvers)
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(α[i][j][range2plot]*180/π, cl[i][j][range2plot], c=colors[i], lw=lw, ls=lstyles[j], label=string(aeroSolver.name," - ",elemLabels[j]))
    end
end
display(plt_cla)
savefig(string(absPath,"/wingDStest_cla.pdf"))

# cm vs α
plt_cma = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legend=:bottomleft, legendfontsize=12)
scatter!(cmRef[1,:], cmRef[2,:], c=:black, ms=ms, label=false)
for (i,aeroSolver) in enumerate(aeroSolvers)
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(α[i][j][range2plot]*180/π, cm[i][j][range2plot], c=colors[i], lw=lw, ls=lstyles[j], label=false)
    end
end
display(plt_cma)
savefig(string(absPath,"/wingDStest_cma.pdf"))

# cd vs α
plt_cda = plot(xlabel="Pitch angle [deg]", ylabel="\$c_d\$", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=12)
scatter!(cdRef[1,:], cdRef[2,:], c=:black, ms=ms, label=false)
for (i,aeroSolver) in enumerate(aeroSolvers)
    ind = findfirst(x-> x>=t[i][end]-lastCycles2plot*τ, t[i])
    range2plot = ind:length(t[i])
    for (j,e) in enumerate(elemRangePlot)
        plot!(α[i][j][range2plot]*180/π, cdrag[i][j][range2plot], c=colors[i], lw=lw, ls=lstyles[j], label=false)
    end
end
display(plt_cda)
savefig(string(absPath,"/wingDStest_cda.pdf"))

println("Finished wingDStestPlotGenerator.jl")