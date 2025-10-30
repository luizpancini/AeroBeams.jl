using Plots, ColorSchemes, JLD2

# Test configurations
tipMassConfig = "LE_m0_o0"
UsweepType = "upsweep"
aeroSolver = Indicial()
tipLossType = "Exponential"
θRange = π/180*vcat(1,3,5,7,10)
URange = [63.5:0.5:74.5,
          50.5:0.5:58.0,
          43.5:0.5:49.5,
          38.5:0.5:43.5,
          34.5:0.5:37.0]

# Set paths
relPath = "/dev/sweptPazy/S0/outputs/PazyWingPitchRangeURangeLCO"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
interactivePlots = false
θcolors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 3
msw = 1
α = 0.5
if interactivePlots
    plotlyjs()
else
    gr()
end          

# --- Plots per root pitch angle for several airspeeds ---
# Sweep root pitch angle
for (i,θ) in enumerate(θRange)
    # --- Initialize plots ---
    # Tip OOP 
    plt_tipOOPAll = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title="\$\\theta=\$ $(round(θ*180/π,digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # Tip AoA
    plt_tipAoAAll = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title="\$\\theta=\$ $(round(θ*180/π,digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # Root LE strains
    plt_rootEpsLEAll = plot(xlabel="Time [s]", ylabel="Root LE strains (\$\\mu\$)", title="\$\\theta=\$ $(round(θ*180/π,digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # Root TE strains
    plt_rootEpsTEAll = plot(xlabel="Time [s]", ylabel="Root TE strains (\$\\mu\$)", title="\$\\theta=\$ $(round(θ*180/π,digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # Set colors
    Ucolors = cgrad(:rainbow, length(URange[i]), categorical=true)
    # Sweep airspeed
    for (j,U) in enumerate(URange[i])
        # Load data
        idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_tipOOP.jld2" tipOOP
        @load absPath*"/"*idString*"_tipAoA.jld2" tipAoA
        @load absPath*"/"*idString*"_tipTwist.jld2" tipTwist
        @load absPath*"/"*idString*"_rootEpsLE.jld2" rootEpsLE
        @load absPath*"/"*idString*"_rootEpsTE.jld2" rootEpsTE
        # Plots
        plot!(plt_tipOOPAll, t, tipOOP/L*100, lw=lw, c=Ucolors, lz=U, label=false)
        plot!(plt_tipAoAAll, t, tipAoA*180/π, lw=lw, c=Ucolors, lz=U, label=false)
        plot!(plt_rootEpsLEAll, t, rootEpsLE*1e6, lw=lw, c=Ucolors, lz=U, label=false)
        plot!(plt_rootEpsTEAll, t, rootEpsTE*1e6, lw=lw, c=Ucolors, lz=U, label=false)
    end
    # Display and save plots
    display(plt_tipOOPAll)
    display(plt_tipAoAAll)
    display(plt_rootEpsLEAll)
    display(plt_rootEpsTEAll)
    savefig(plt_tipOOPAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipOOP_all.pdf"))
    savefig(plt_tipAoAAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipAoA_all.pdf"))
    savefig(plt_rootEpsLEAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_rootEpsLE_all.pdf"))
    savefig(plt_rootEpsTEAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_rootEpstE_all.pdf"))
end

# --- Bifurcation diagrams ---
# Time span for bifurcation diagrams
time_BD = 1

# Bifurcation diagrams: tip OOP
plt_tipOOP_BD = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[30,80], ylims=[10,30], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for (j,U) in enumerate(URange[i])
        # Case ID
        idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        # Load data
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_tipOOP.jld2" tipOOP
        # Compute extrema over last time_BD seconds
        ind = findfirst(x -> x >= t[end]-time_BD, t)
        ext = local_extrema(tipOOP[ind:end]*100/L)
        # Plot
        scatter!(U*ones(length(ext.maxVal)), ext.maxVal, c=θcolors[i], shape=:utriangle, ms=ms, msw=0, label=false)
        scatter!(U*ones(length(ext.minVal)), ext.minVal, c=θcolors[i], shape=:dtriangle, ms=ms, msw=0, label=false)
    end
end
display(plt_tipOOP_BD)
savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipOOP_BD.pdf"))

# Bifurcation diagrams: tip AoA
plt_tipAoA_BD = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[30,80], ylims=[-5,15], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for (j,U) in enumerate(URange[i])
        # Case ID
        idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        # Load data
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_tipAoA.jld2" tipAoA
        # Compute extrema over last time_BD seconds
        ind = findfirst(x -> x >= t[end]-time_BD, t)
        ext = local_extrema(tipAoA[ind:end]*180/π)
        # Plot
        scatter!(U*ones(length(ext.maxVal)), ext.maxVal, c=θcolors[i], shape=:utriangle, ms=ms, msw=0, label=false)
        scatter!(U*ones(length(ext.minVal)), ext.minVal, c=θcolors[i], shape=:dtriangle, ms=ms, msw=0, label=false)
    end
end
display(plt_tipAoA_BD)
savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipAoA_BD.pdf"))

# Bifurcation diagrams: root LE strains
plt_rootEpsLE_BD = plot(xlabel="Airspeed [m/s]", ylabel="Root LE strains (\$\\mu\$)", xlims=[30,80], ylims=[0,3000], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for (j,U) in enumerate(URange[i])
        # Case ID
        idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        # Load data
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_rootEpsLE.jld2" rootEpsLE
        # Compute extrema over last time_BD seconds
        ind = findfirst(x -> x >= t[end]-time_BD, t)
        ext = local_extrema(rootEpsLE[ind:end]*1e6)
        # Plot
        scatter!(U*ones(length(ext.maxVal)), ext.maxVal, c=θcolors[i], shape=:utriangle, ms=ms, msw=0, label=false)
        scatter!(U*ones(length(ext.minVal)), ext.minVal, c=θcolors[i], shape=:dtriangle, ms=ms, msw=0, label=false)
    end
end
display(plt_rootEpsLE_BD)
savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_rootEpsLE_BD.pdf"))

# Bifurcation diagrams: root TE strains
plt_rootEpsTE_BD = plot(xlabel="Airspeed [m/s]", ylabel="Root TE strains (\$\\mu\$)", xlims=[30,80], ylims=[0,3000], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for (j,U) in enumerate(URange[i])
        # Case ID
        idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        # Load data
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_rootEpsTE.jld2" rootEpsTE
        # Compute extrema over last time_BD seconds
        ind = findfirst(x -> x >= t[end]-time_BD, t)
        ext = local_extrema(rootEpsTE[ind:end]*1e6)
        # Plot
        scatter!(U*ones(length(ext.maxVal)), ext.maxVal, c=θcolors[i], shape=:utriangle, ms=ms, msw=0, label=false)
        scatter!(U*ones(length(ext.minVal)), ext.minVal, c=θcolors[i], shape=:dtriangle, ms=ms, msw=0, label=false)
    end
end
display(plt_rootEpsTE_BD)
savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_rootEpsTE_BD.pdf"))

println("Finished PazyWingPitchRangeURangeLCOPlotGenerator.jl")