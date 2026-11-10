using AeroBeams, JLD2

# Test configurations
tipMassConfig = "LE_m0_o0"
UsweepType = "upsweep"
aeroSolver = Indicial()
tipLossType = "Exponential"
θRange = deg2rad.([1,3,5,7,10])
URange = [63.5:0.5:74.5,
          50.5:0.5:58.0,
          43.5:0.5:49.5,
          38.5:0.5:43.5,
          34.5:0.5:37.0]

# Fixed geometrical properties
_,L = geometrical_properties_Pazy()

## Plots
# Set paths
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingPitchRangeURangeLCO"
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
    AeroBeams.Plots.plotlyjs()
else
    gr()
end

# --- Plots per root pitch angle for several airspeeds ---
# Sweep root pitch angle
for (i,θ) in enumerate(θRange)
    # --- Initialize plots ---
    # Tip OOP
    plt_tipOOPAll = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title="\$\\theta=\$ $(round(rad2deg(θ),digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # Tip AoA
    plt_tipAoAAll = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title="\$\\theta=\$ $(round(rad2deg(θ),digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # 1/6-span strains
    plt_epsAll = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\$)", title="\$\\theta=\$ $(round(rad2deg(θ),digits=1)) deg", tickfont=font(ts), guidefont=font(fs), colorbar_title="Airspeed [m/s]")
    # Set colors
    Ucolors = cgrad(:rainbow, length(URange[i]), categorical=true)
    # Sweep airspeed
    for (j,U) in enumerate(URange[i])
        # Load data
        idString = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_tipOOP.jld2" tipOOP
        @load absPath*"/"*idString*"_tipAoA.jld2" tipAoA
        @load absPath*"/"*idString*"_tipTwist.jld2" tipTwist
        @load absPath*"/"*idString*"_epsSG.jld2" epsSG
        # Plots
        plot!(plt_tipOOPAll, t, tipOOP/L*100, lw=lw, c=Ucolors, lz=U, label=false)
        plot!(plt_tipAoAAll, t, rad2deg.(tipAoA), lw=lw, c=Ucolors, lz=U, label=false)
        plot!(plt_epsAll, t, epsSG*1e6, lw=lw, c=Ucolors, lz=U, label=false)
    end
    # Display and save plots
    display(plt_tipOOPAll)
    display(plt_tipAoAAll)
    display(plt_epsAll)
    savefig(plt_tipOOPAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipOOP_all.pdf"))
    savefig(plt_tipAoAAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipAoA_all.pdf"))
    savefig(plt_epsAll,string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_eps_all.pdf"))
end

# --- Bifurcation diagrams ---
# Time span for bifurcation diagrams
time_BD = 1

# Bifurcation diagrams: tip OOP
plt_tipOOP_BD = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[30,80], ylims=[10,30], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for (j,U) in enumerate(URange[i])
        # Case ID
        idString = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
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
        idString = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        # Load data
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_tipAoA.jld2" tipAoA
        # Compute extrema over last time_BD seconds
        ind = findfirst(x -> x >= t[end]-time_BD, t)
        ext = local_extrema(rad2deg.(tipAoA[ind:end]))
        # Plot
        scatter!(U*ones(length(ext.maxVal)), ext.maxVal, c=θcolors[i], shape=:utriangle, ms=ms, msw=0, label=false)
        scatter!(U*ones(length(ext.minVal)), ext.minVal, c=θcolors[i], shape=:dtriangle, ms=ms, msw=0, label=false)
    end
end
display(plt_tipAoA_BD)
savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_tipAoA_BD.pdf"))

# Bifurcation diagrams: 1/6-span strains
plt_eps_BD = plot(xlabel="Airspeed [m/s]", ylabel="1/6-span strains (\$\\mu\$)", xlims=[30,80], ylims=[0,3000], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    for (j,U) in enumerate(URange[i])
        # Case ID
        idString = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name, "_", UsweepType, "_U", round(U,digits=3))
        # Load data
        @load absPath*"/"*idString*"_t.jld2" t
        @load absPath*"/"*idString*"_epsSG.jld2" epsSG
        # Compute extrema over last time_BD seconds
        ind = findfirst(x -> x >= t[end]-time_BD, t)
        ext = local_extrema(epsSG[ind:end]*1e6)
        # Plot
        scatter!(U*ones(length(ext.maxVal)), ext.maxVal, c=θcolors[i], shape=:utriangle, ms=ms, msw=0, label=false)
        scatter!(U*ones(length(ext.minVal)), ext.minVal, c=θcolors[i], shape=:dtriangle, ms=ms, msw=0, label=false)
    end
end
display(plt_eps_BD)
savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",tipLossType, "_", aeroSolver.name,"_eps_BD.pdf"))

println("Finished PazyWingPitchRangeURangeLCOPlotGenerator.jl")
