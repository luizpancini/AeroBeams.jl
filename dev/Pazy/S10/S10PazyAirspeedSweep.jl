using AeroBeams

# Root pitch angle offset (assumed) [rad]
θoffset = deg2rad(0.0)

# Root pitch angle [rad]
θ = deg2rad(10) + θoffset

# Profile for airspeed sweep [m/s]
τ = 30
Umin = 40
Umax = 44
U = t -> ifelse(t<=τ, Umin + (Umax-Umin)*t/τ, Umax)

# Flag to solve flutter problem
solveFlutter = true

# Airspeed range for flutter problem [m/s]
URange = collect(0.5:0.5:90)

# Sweep angle [rad]
Λ = deg2rad(10)

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag to update tip correction with airspeed
tipLossFunctionIsAirspeedDependent = true

# Aerodynamic solver
aeroSolver = BLi()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(Λ=Λ)

# Tip mass configuration
tipMassConfig = :LE

# Time variables
Δt = 1e-4
tf = τ
trackingFrequency = clamp(round(Int,1e-3/Δt), 1, 10)

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 40

# Set system solver options
systemSolver = aeroSolver isa DynamicStallSolver ?
               create_NewtonRaphson(displayStatus=false) :
               create_NLsolver(method=:newton,linesearch=false,displayStatus=false)

# Strain gauge coordinates (the spar cs is 60 x 2.25 mm)
SGelem = 3          # element number where the strain gauge is located (see Fig. 12 of https://doi.org/10.2514/1.J060621: FBG^1_3 is @ 127 mm spanwise, but wing root is at 32 mm -> actual spanwise position is 127 - 32 = 95 mm which is at the midpoint of element 3)
x2SG = 50e-3/2      # chordwise position: close to LE
x3SG = 2.25e-3/2    # height: on top

# Case ID
idString = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_Umax", Umax, "_", tipLossType, "_", aeroSolver.name)

# Set paths
relPath = "/dev/Pazy/S10/outputs/S10PazyAirspeedSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

## --- Solve flutter problem ---
if solveFlutter

    println("Solving flutter problem")

    # Initialize outputs
    flutterProblem = Array{EigenProblem}(undef,length(URange))
    eigensols = Array{EigenSolution}(undef,length(URange))
    tipOOP_f = Array{Float64}(undef,length(URange))
    tipAoA_f = Array{Float64}(undef,length(URange))
    epsSG_f = Array{Float64}(undef,length(URange))

    # Initialize model
    flutterModel,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)

    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Update airspeed
        update_v_A!(flutterModel; v_A=[0,U,0])
        # Set initial guess solution as previously available one
        x0 = j == 1 ? zeros(0) : flutterProblem[j-1].systemArrays.x
        # Create and solve problem
        flutterProblem[j] = create_EigenProblem(model=flutterModel,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(flutterProblem[j])
        # Eigensolution
        eigensols[j] = flutterProblem[j].eigensolution
        # Kinematics
        states = flutterProblem[j].loadstepHistory.states
        aeroVars = flutterProblem[j].loadstepHistory.aeroVariables
        tipOOP_f[j] = -states[end][nElem].nodal.u2[1]
        tipAoA_f[j] = aeroVars[end][nElem].unsteadyKin.αₑ
        ϵ11SG_f = states[end][SGelem].midpoint.γ[1]
        κ2SG_f = states[end][SGelem].midpoint.κ[2]
        κ3SG_f = states[end][SGelem].midpoint.κ[3]
        epsSG_f[j] = (ϵ11SG_f .- κ2SG_f*x3SG .- κ3SG_f*x2SG)
    end

    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols)

    # Get frequencies and dampings of oscillatory modes
    freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1), fmax=Hz2rads(50), maxTracks=nModes, dropOutsideBounds=true)

    # Damping ratios
    dampingRatios = dampsOsc./freqsOsc

    # Flutter onset and offset speeds, frequencies and kinematics
    flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetModes,flutterOnsetTipOOPs,flutterOnsetTipAoAs,flutterOnsetEpss,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetModes,flutterOffsetTipOOPs,flutterOffsetTipAoAs,flutterOffsetEpss = find_flutter(URange,freqsOsc,dampsOsc,tipOOP_f,tipAoA_f,epsSG_f)
    if !isempty(flutterOnsetSpeeds)
        flutterOnsetSpeed = flutterOnsetSpeeds[1]
        flutterOnsetFreq = flutterOnsetFreqs[1]
        flutterOnsetTipOOP = flutterOnsetTipOOPs[1]
        flutterOnsetTipAoA = flutterOnsetTipAoAs[1]
        flutterOnsetEps = flutterOnsetEpss[1]
        println("Flutter onset speed = $(round(flutterOnsetSpeed,digits=1)) m/s, frequency = $(round(rads2Hz(flutterOnsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOnsetTipAoA),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOnsetEps*1e6)) microns")
    else
        println("Flutter not found")
    end
    if !isempty(flutterOffsetSpeeds)
        flutterOffsetSpeed = flutterOffsetSpeeds[1]
        flutterOffsetFreq = flutterOffsetFreqs[1]
        flutterOffsetTipOOP = flutterOffsetTipOOPs[1]
        flutterOffsetTipAoA = flutterOffsetTipAoAs[1]
        flutterOffsetEps = flutterOffsetEpss[1]
        println("Flutter offset speed = $(round(flutterOffsetSpeed,digits=1)) m/s, frequency = $(round(rads2Hz(flutterOffsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOffsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOffsetTipAoA),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOffsetEps*1e6)) microns")
    end

    ## Plots
    modeColors = cgrad(:rainbow, nModes, categorical=true)
    ts = 10
    fs = 14
    lw = 2
    gr()
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    for mode in 1:nModes
        plot!(URange, rads2Hz.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URange, dampingRatios[mode,:], c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(plt_Vgf,string(absPath,"/S10PazyAirspeedSweep_",idString,"_Vgf.pdf"))
end

## --- Solve dynamic problem ---
println("Solving dynamic problem")

# Model
dynModel,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,airspeed=U,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent,sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)

# Create and solve steady problem for initial airspeed
steadyProblem = create_SteadyProblem(model=dynModel)
solve!(steadyProblem)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=dynModel,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x,trackingFrequency=trackingFrequency)
solve!(dynamicProblem)

# Unpack numerical solution
sol = dynamicProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
tipAoA = [aeroVars[i][nElem].unsteadyKin.αₑ for i in eachindex(t)]
tipOOP = -[states[i][nElem].nodal.u2[1] for i in eachindex(t)]
ϵ11SG = [states[i][SGelem].midpoint.γ[1] for i in eachindex(t)]
κ2SG = [states[i][SGelem].midpoint.κ[2] for i in eachindex(t)]
κ3SG = [states[i][SGelem].midpoint.κ[3] for i in eachindex(t)]
epsSG = (ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG)

## Plots
# Animation
plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=round(Int,5e-2/Δt),plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/S10PazyAirspeedSweep_",idString,"_deformation.gif"),displayProgress=true)

# Plot configurations
showAirspeedInPlots = true
interactivePlots = false
colors = cgrad(:rainbow, 2, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 6
msw = 1
if interactivePlots
    plotlyjs()
else
    gr()
end

# Airspeed
if solveFlutter
    jF = findfirst(ti -> U(ti) > flutterOnsetSpeed, t)
end
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, U.(t), lw=lw, c=:black, label=false)
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], [U(t[jF])], marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], [Umin, U(t[jF])], ls=:dash, c=:red, lw=lw, label=false)
    plot!([0, t[jF]], [U(t[jF]), U(t[jF])], ls=:dash, c=:red, lw=lw, label=false)
end
display(plt_U)
savefig(string(absPath,"/S10PazyAirspeedSweep_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, tipOOP/L*100, lw=lw, c=colors[1], label="AeroBeams")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=colors[2], label="Revivo & Raveh (2025)")
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], [tipOOP[jF]]/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], [tipOOP[1], tipOOP[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, t[jF]], [tipOOP[jF], tipOOP[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots
    plt_tipOOP_twin = AeroBeams.Plots.twinx(plt_tipOOP)
    plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), extrema(t))
    plot!(plt_tipOOP_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipOOP)
savefig(string(absPath,"/S10PazyAirspeedSweep_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, rad2deg.(tipAoA), lw=lw, c=colors[1], label=false)
plot!([NaN], [NaN], lw=lw, ls=:dash, c=colors[2], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], rad2deg.([tipAoA[jF]]), marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], rad2deg.([tipAoA[1], tipAoA[jF]]), ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, t[jF]], rad2deg.([tipAoA[jF], tipAoA[jF]]), ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots
    plt_tipAOA_twin = AeroBeams.Plots.twinx(plt_tipAOA)
    plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_tipAOA_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/S10PazyAirspeedSweep_",idString,"_tipAoA.pdf"))

# 1/6-span strains
plt_eps = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\\epsilon\$)", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, epsSG*1e6, lw=lw, c=colors[1], alpha=0.7, label=false)
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], [epsSG[jF]]*1e6, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], [epsSG[1], epsSG[jF]]*1e6, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, t[jF]], [epsSG[jF], epsSG[jF]]*1e6, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots
    plt_eps_twin = AeroBeams.Plots.twinx(plt_eps)
    plot!(plt_eps_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_eps_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_eps)
savefig(string(absPath,"/S10PazyAirspeedSweep_",idString,"_eps.pdf"))

# 1/6-span strains spectrogram
tBegin = 0
jBegin = findfirst(x-> x>=tBegin, t)
tCropped = t[jBegin:end] .- t[jBegin]
epsSGCropped = epsSG[jBegin:end]
f_eps, t_eps, S_eps = AeroBeams.spectrogram(tCropped, epsSGCropped*1e6, window=2^12)
U_vec = U.(tCropped)
U_eps = U_vec[searchsortedfirst.(Ref(tCropped), t_eps)]
plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title="1/6-span strains spectrogram", xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
plot!([U_vec[1], U_vec[end]], rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
plot!([U_vec[1], U_vec[end]], 2*rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
plot!([U_vec[1], U_vec[end]], 3*rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
display(plt_eps_spectrum)
savefig(string(absPath,"/S10PazyAirspeedSweep_",idString,"_epsSpectrogram.pdf"))

println("Finished S10PazyAirspeedSweep.jl")
