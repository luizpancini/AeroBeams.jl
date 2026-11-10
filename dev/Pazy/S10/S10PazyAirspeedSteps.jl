using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(10)

# Root pitch angle offset (assumed) [rad]
θoffset = deg2rad(0)

# Root pitch angle [rad]
θ = deg2rad(10) + θoffset

# Flag to solve flutter problem
solveFlutter = true

# Airspeed range for flutter problem [m/s]
URangeFlutter = collect(0.5:0.5:90)

# Flag to model the wingtip rod as flexible (true) or rigid (false)
flexibleWingtipRod = true

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = flexibleWingtipRod

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

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 40

# Strain gauge coordinates (the spar cs is 60 x 2.25 mm)
SGelem = 3          # element number where the strain gauge is located (see Fig. 12 of https://doi.org/10.2514/1.J060621: FBG^1_3 is @ 127 mm spanwise, but wing root is at 32 mm -> actual spanwise position is 127 - 32 = 95 mm which is at the midpoint of element 3)
x2SG = 50e-3/2      # chordwise position: close to LE
x3SG = 2.25e-3/2    # height: on top

# Set system solver options
systemSolver = aeroSolver isa DynamicStallSolver ?
               create_NewtonRaphson(absTol=1e-7,displayStatus=false) :
               create_NLsolver(method=:newton,linesearch=false,displayStatus=false)

# Case ID
idString = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name)

# Set paths
relPath = "/dev/Pazy/S10/outputs/S10PazyAirspeedSteps"
absPath = string(pwd(),relPath)
mkpath(absPath)

## --- Solve flutter problem ---
if solveFlutter
    println("Solving flutter problem")

    # Initialize outputs
    flutterProblem = Array{EigenProblem}(undef,length(URangeFlutter))
    eigensols = Array{EigenSolution}(undef,length(URangeFlutter))
    tipOOP_f = Array{Float64}(undef,length(URangeFlutter))
    tipAoA_f = Array{Float64}(undef,length(URangeFlutter))
    epsSG_f = Array{Float64}(undef,length(URangeFlutter))

    # Initialize model
    flutterModel,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)

    # Sweep airspeed
    for (j,U) in enumerate(URangeFlutter)
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
    flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetModes,flutterOnsetTipOOPs,flutterOnsetTipAoAs,flutterOnsetEpss,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetModes,flutterOffsetTipOOPs,flutterOffsetTipAoAs,flutterOffsetEpss = find_flutter(URangeFlutter,freqsOsc,dampsOsc,tipOOP_f,tipAoA_f,epsSG_f)
    if !isempty(flutterOnsetSpeeds)
        flutterOnsetSpeed = flutterOnsetSpeeds[1]
        flutterOnsetFreq = flutterOnsetFreqs[1]
        flutterOnsetTipOOP = flutterOnsetTipOOPs[1]
        flutterOnsetTipAoA = flutterOnsetTipAoAs[1]
        flutterOnsetEps = flutterOnsetEpss[1]
        println("Flutter onset speed = $(round(flutterOnsetSpeed,digits=2)) m/s, frequency = $(round(rads2Hz(flutterOnsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOnsetTipAoA),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOnsetEps*1e6)) microns")
    else
        println("Flutter not found")
    end
    if !isempty(flutterOffsetSpeeds)
        flutterOffsetSpeed = flutterOffsetSpeeds[1]
        flutterOffsetFreq = flutterOffsetFreqs[1]
        flutterOffsetTipOOP = flutterOffsetTipOOPs[1]
        flutterOffsetTipAoA = flutterOffsetTipAoAs[1]
        flutterOffsetEps = flutterOffsetEpss[1]
        println("Flutter offset speed = $(round(flutterOffsetSpeed,digits=2)) m/s, frequency = $(round(rads2Hz(flutterOffsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOffsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOffsetTipAoA),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOffsetEps*1e6)) microns")
    end

    ## Plots
    modeColors = cgrad(:rainbow, nModes, categorical=true)
    ts = 10
    fs = 16
    lw = 2
    gr()
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URangeFlutter), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    for mode in 1:nModes
        plot!(URangeFlutter, rads2Hz.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URangeFlutter), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    plot!(URangeFlutter, zeros(length(URangeFlutter)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URangeFlutter, dampingRatios[mode,:], c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(plt_Vgf,string(absPath,"/S10PazyAirspeedSteps_",idString,"_Vgf.pdf"))
end

## --- Dynamic problems setup ---

# Airspeed range for dynamic problems [m/s]
ΔU = 0.25
# URangeDyn = collect(floor((flutterOnsetSpeed-ΔU)/ΔU)*ΔU:ΔU:ceil((flutterOffsetSpeed+ΔU)/ΔU)*ΔU)
URangeDyn = collect(floor((flutterOnsetSpeed-ΔU)/ΔU)*ΔU:ΔU:ceil((flutterOnsetSpeed+2*ΔU)/ΔU)*ΔU)

# Airspeed profile variables
N = length(URangeDyn)
τ₀Range = 1*ones(N-1)  # Steady lower value time
τ₁Range = 1*ones(N-1)  # Linear transition time
τ₂Range = 2*ones(N-1)  # Steady upper value time
# τ₂Range[findfirst(x->x>=flutterOnsetSpeed, URangeDyn)] = 10

# Time variables
Δt = 1e-4
trackingFrequency = clamp(round(Int,1e-3/Δt), 1, 10)

# --- Solve initial steady problem ---
println("Solving initial steady problem")

# Create and solve steady problem for initial airspeed
steadyModel,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,airspeed=URangeDyn[1],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent,sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)
steadyProblem = create_SteadyProblem(model=steadyModel)
solve!(steadyProblem)

# --- Solve dynamic problems ---
println("Solving dynamic problems")

# Initialize outputs
dynamicProblem = Array{DynamicProblem}(undef,N-1)
Uvec = Array{Vector{Float64}}(undef,N-1)
t = Array{Vector{Float64}}(undef,N-1)
tipOOP = Array{Vector{Float64}}(undef,N-1)
tipAoA = Array{Vector{Float64}}(undef,N-1)
epsSG = Array{Vector{Float64}}(undef,N-1)
timespan = Array{Float64}(undef,N-1)
initialTime = Array{Float64}(undef,N-1)
finalTime = Array{Float64}(undef,N-1)

# Step over airspeeds
for (j,U) in enumerate(URangeDyn[1:end-1])
    # Display progress
    println("Solving for U = $U -> $(URangeDyn[j+1]) m/s")
    # Set initial solution
    x0 = j == 1 ? steadyProblem.systemArrays.x : dynamicProblem[j-1].systemArrays.x
    # Span, initial and final times of simulation
    timespan[j] = τ₀Range[j] + τ₁Range[j] + τ₂Range[j]
    initialTime[j] = j == 1 ? 0.0 : finalTime[j-1]
    finalTime[j] = initialTime[j] + timespan[j]
    # Set airspeed profile
    airspeed = t -> ifelse(t<=initialTime[j]+τ₀Range[j], U, ifelse(t>=initialTime[j]+τ₀Range[j]+τ₁Range[j], URangeDyn[j+1], U + (URangeDyn[j+1]-U)*(t-τ₀Range[j]-initialTime[j])/τ₁Range[j]))
    # Model
    dynModel,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,pitchAngle=θ,airspeed=airspeed,tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent,sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)
    # Reference dynamic problem
    refDynamicProblem = j==1 ? nothing : dynamicProblem[j-1]
    # Create and solve dynamic problem
    dynamicProblem[j] = create_DynamicProblem(model=dynModel,initialTime=initialTime[j],finalTime=finalTime[j],Δt=Δt,systemSolver=systemSolver,x0=x0,refDynamicProblem=refDynamicProblem,trackingFrequency=trackingFrequency)
    solve!(dynamicProblem[j])
    # Unpack numerical solution
    sol = dynamicProblem[j].solutionHistory
    # Skip the initial sample of continuation segments to avoid duplicated time points
    idx = j == 1 ? eachindex(sol.t) : 2:length(sol.t)
    t[j] = sol.t[idx]
    Uvec[j] = airspeed.(t[j])
    tipAoA[j] = [sol.aeroVariables[i][nElem].unsteadyKin.αₑ for i in idx]
    tipOOP[j] = -[sol.states[i][nElem].nodal.u2[1] for i in idx]
    ϵ11SG = [sol.states[i][SGelem].midpoint.γ[1] for i in idx]
    κ2SG = [sol.states[i][SGelem].midpoint.κ[2] for i in idx]
    κ3SG = [sol.states[i][SGelem].midpoint.κ[3] for i in idx]
    epsSG[j] = (ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG)
end

# Concatenate dynamic problems outputs
tAll = vcat(t...)
UAll = vcat(Uvec...)
tipAoAAll = vcat(tipAoA...)
tipOOPAll = vcat(tipOOP...)
epsSGAll = vcat(epsSG...)

# Set case ID for dynamic plots
idString = URangeDyn[end] > URangeDyn[1] ? string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name, "_upsteps") : string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name, "_downsteps")

## Plots
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
α = 0.5
if interactivePlots
    plotlyjs()
else
    gr()
end

# Airspeed
if solveFlutter
    jF = findfirst(x -> x >= flutterOnsetSpeed, UAll)
end
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", xlims=extrema(tAll), ylims=extrema(UAll), tickfont=font(ts), guidefont=font(fs))
plot!(tAll, UAll, lw=lw, c=:black, label=false)
if solveFlutter && !isnothing(jF)
    scatter!([tAll[jF]], [UAll[jF]], marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([tAll[jF], tAll[jF]], [URangeDyn[1], UAll[jF]], ls=:dash, c=:red, lw=lw, label=false)
    plot!([0, tAll[jF]], [UAll[jF], UAll[jF]], ls=:dash, c=:red, lw=lw, label=false)
end
display(plt_U)
savefig(string(absPath,"/S10PazyAirspeedSteps_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=extrema(tAll), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(tAll, tipOOPAll/L*100, lw=lw, c=colors[1], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([tAll[jF]], [tipOOPAll[jF]]/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([tAll[jF], tAll[jF]], [tipOOPAll[1], tipOOPAll[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, tAll[jF]], [tipOOPAll[jF], tipOOPAll[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_tipOOP_twin = AeroBeams.Plots.twinx(plt_tipOOP)
    plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(tAll), ylims=extrema(UAll))
    plot!(plt_tipOOP_twin, tAll, UAll, lw=lw, c=:black, label=false)
end
display(plt_tipOOP)
savefig(string(absPath,"/S10PazyAirspeedSteps_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=extrema(tAll), tickfont=font(ts), guidefont=font(fs))
plot!(tAll, rad2deg.(tipAoAAll), lw=lw, c=colors[1], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([tAll[jF]], rad2deg.([tipAoAAll[jF]]), marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([tAll[jF], tAll[jF]], rad2deg.([tipAoAAll[1], tipAoAAll[jF]]), ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, tAll[jF]], rad2deg.([tipAoAAll[jF], tipAoAAll[jF]]), ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_tipAOA_twin = AeroBeams.Plots.twinx(plt_tipAOA)
    plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(tAll), ylims=extrema(UAll))
    plot!(plt_tipAOA_twin, tAll, UAll, lw=lw, c=:black, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/S10PazyAirspeedSteps_",idString,"_tipAoA.pdf"))

# 1/6-span strains
plt_eps = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\\epsilon\$)", xlims=extrema(tAll), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(tAll, epsSGAll*1e6, lw=lw, c=colors[1], alpha=α, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_eps_twin = AeroBeams.Plots.twinx(plt_eps)
    plot!(plt_eps_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(tAll), ylims=extrema(UAll))
    plot!(plt_eps_twin, tAll, UAll, lw=lw, c=:black, label=false)
end
display(plt_eps)
savefig(string(absPath,"/S10PazyAirspeedSteps_",idString,"_eps.pdf"))

# 1/6-span strains spectrogram
tBegin = 0
jBegin = findfirst(x-> x>=tBegin, tAll)
tCropped = tAll[jBegin:end] .- tAll[jBegin]
epsSGCropped = epsSGAll[jBegin:end]
f_eps, t_eps, S_eps = AeroBeams.spectrogram(tCropped, epsSGCropped*1e6, window=2^12)
U_eps = UAll[searchsortedfirst.(Ref(tAll), t_eps)]
plt_eps_spectrum = heatmap(U_eps, f_eps, S_eps; ylims=[0,100], title="1/6-span strains spectrogram", xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
plot!([flutterOnsetSpeed, flutterOnsetSpeed], [0, 100], ls=:dash, c=:black, lw=lw, label=false)
plot!([UAll[1], U_eps[end]], rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
plot!([UAll[1], U_eps[end]], 2*rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
plot!([UAll[1], U_eps[end]], 3*rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
display(plt_eps_spectrum)
savefig(string(absPath,"/S10PazyAirspeedSteps_",idString,"_epsSpectrogram.pdf"))

println("Finished S10PazyAirspeedSteps.jl")
