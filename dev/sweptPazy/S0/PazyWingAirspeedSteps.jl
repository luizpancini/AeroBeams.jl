using AeroBeams, LinearAlgebra, DelimitedFiles, LinearInterpolations, Plots, ColorSchemes

# Flag to solve flutter problem
solveFlutter = true

# Select tip mass configuration (set as XX_mY_oZ, where XX is either LE or TE, Y is the mass in grams and Z is the offset in cm)
tipMassConfig = "LE_m0_o0"

# Root pitch angle offset (assumed)
θoffset = 0.0*π/180

# Root pitch angle
θ = 7π/180 + θoffset

# Airspeed range for flutter problem
URangeFlutter = collect(20:0.25:90)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "Exponential"

# Flag to update tip correction with airspeed
tipLossFunctionIsAirspeedDependent = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Set tip mass configuration variables
m = match(r"^(LE|TE)_m(\d+(?:\.\d+)?)_o(\d+(?:\.\d+)?)$", tipMassConfig)
@assert m !== nothing "Invalid tipMassConfig string"
pos = m.captures[1]
tipMass = 1e-3*parse(Float64, m.captures[2])
tipMassPosOffset = 1e-2*parse(Float64, m.captures[3])
tipMassPos = pos == "LE" ? chord*normSparPos + tipMassPosOffset : -(chord*(1-normSparPos) + tipMassPosOffset)

# Number of modes for flutter computation
nModes = 3

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
ySG_LE = 50e-3/2    # y-position: on LE
ySG_TE = -50e-3/2   # y-position: on TE
zSG = 2.25e-3/2     # z-position: on top

# Case ID
idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name)

# --- Solve flutter problem ---
if solveFlutter
    println("Solving flutter problem")

    # Initialize outputs
    flutterProblem = Array{EigenProblem}(undef,length(URangeFlutter))
    untrackedFreqs = Array{Vector{Float64}}(undef,length(URangeFlutter))
    untrackedDamps = Array{Vector{Float64}}(undef,length(URangeFlutter))
    untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URangeFlutter))
    tipOOP_f = Array{Float64}(undef,length(URangeFlutter))
    tipAoA_f = Array{Float64}(undef,length(URangeFlutter))
    rootEps_f = Array{Float64}(undef,length(URangeFlutter))

    # Sweep airspeed
    for (j,U) in enumerate(URangeFlutter)
        # Model
        flutterModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,θ=θ,airspeed=U,tipMass=tipMass,ηtipMass=[0;tipMassPos;0]))
        # Set initial guess solution as previously available one
        x0 = j == 1 ? zeros(0) : flutterProblem[j-1].x
        # Create and solve problem
        flutterProblem[j] = create_EigenProblem(model=flutterModel,nModes=nModes,frequencyFilterLimits=[2π,Inf],x0=x0)
        solve!(flutterProblem[j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[j] = flutterProblem[j].frequenciesOscillatory
        untrackedDamps[j] = flutterProblem[j].dampingsOscillatory
        untrackedEigenvectors[j] = flutterProblem[j].eigenvectorsOscillatoryCplx
        # Kinematics
        tipOOP_f[j] = -flutterProblem[j].nodalStatesOverσ[end][nElem].u_n2[1]
        tipAoA_f[j] = flutterProblem[j].aeroVariablesOverσ[end][nElem].flowAnglesAndRates.αₑ
        ϵ11Root_f = flutterProblem[j].compElementalStatesOverσ[end][1].γ[1]
        κ2Root_f = flutterProblem[j].compElementalStatesOverσ[end][1].κ[2]
        rootEps_f[j] = (ϵ11Root_f .- κ2Root_f*zSG)
    end

    # Apply mode tracking
    freqs,damps,_ = mode_tracking_hungarian(URangeFlutter,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

    # Separate frequencies and damping ratios by mode
    modeFrequencies = Array{Vector{Float64}}(undef,nModes)
    modeDampings = Array{Vector{Float64}}(undef,nModes)
    modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
    for mode in 1:nModes
        modeFrequencies[mode] = [freqs[j][mode] for j in eachindex(URangeFlutter)]
        modeDampings[mode] = [damps[j][mode] for j in eachindex(URangeFlutter)]
        modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
    end

    # Compute flutter variables
    flutterOnsetSpeedOfMode = [fill(0.0,0) for _ in 1:nModes]
    flutterOffsetSpeedOfMode = [fill(0.0,0) for _ in 1:nModes]
    flutterOnsetFreqOfMode = [fill(0.0,0) for _ in 1:nModes]
    flutterOnsetTipOOPOfMode = [fill(0.0,0) for _ in 1:nModes]
    flutterOnsetTipAoAOfMode = [fill(0.0,0) for _ in 1:nModes]
    flutterOnsetRootEpsOfMode = [fill(0.0,0) for _ in 1:nModes]
    for mode in 1:nModes
        jOnset = findall(j -> modeDampings[mode][j] < 0 && modeDampings[mode][j+1] > 0, 1:length(URangeFlutter)-1)
        jOffset = findall(j -> modeDampings[mode][j] > 0 && modeDampings[mode][j+1] < 0, 1:length(URangeFlutter)-1)
        if !isempty(jOnset)
            for jO in jOnset
                push!(flutterOnsetSpeedOfMode[mode],interpolate(modeDampings[mode][jO:jO+1],URangeFlutter[jO:jO+1], 0))
                push!(flutterOnsetFreqOfMode[mode],interpolate(modeDampings[mode][jO:jO+1],modeFrequencies[mode][jO:jO+1], 0))
                push!(flutterOnsetTipOOPOfMode[mode], tipOOP_f[jO])
                push!(flutterOnsetTipAoAOfMode[mode], tipAoA_f[jO])
                push!(flutterOnsetRootEpsOfMode[mode], rootEps_f[jO])
            end
        end
        if !(isempty(jOffset) || isempty(jOnset))
            for jO in jOffset
                push!(flutterOffsetSpeedOfMode[mode],interpolate(-modeDampings[mode][jO:jO+1],URangeFlutter[jO:jO+1],0))
            end
        end
    end
    flutterOnsetSpeedsAll = vcat(filter(!isempty,flutterOnsetSpeedOfMode)...)
    flutterOffsetSpeedsAll = vcat(filter(!isempty,flutterOffsetSpeedOfMode)...)
    flutterOnsetFreqAll = vcat(filter(!isempty,flutterOnsetFreqOfMode)...)
    flutterOnsetTipOOPAll = vcat(filter(!isempty,flutterOnsetTipOOPOfMode)...)
    flutterOnsetTipAoAAll = vcat(filter(!isempty,flutterOnsetTipAoAOfMode)...)
    flutterOnsetRootEpsAll = vcat(filter(!isempty,flutterOnsetRootEpsOfMode)...)
    if !isempty(flutterOnsetSpeedsAll)
        iOnLowest = sortperm(flutterOnsetSpeedsAll)[1]
        flutterOnsetSpeed = flutterOnsetSpeedsAll[iOnLowest]
        flutterOnsetFreq = flutterOnsetFreqAll[iOnLowest]
        flutterOnsetTipOOP = flutterOnsetTipOOPAll[iOnLowest]
        flutterOnsetTipAoA = flutterOnsetTipAoAAll[iOnLowest]
        flutterOnsetRootEps = flutterOnsetRootEpsAll[iOnLowest]
        println("Flutter onset speed = $(round(flutterOnsetSpeed,digits=1)) m/s, frequency = $(round(flutterOnsetFreq/(2*π),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(flutterOnsetTipAoA*180/π,digits=1)) deg, root strains = $(round(Int,flutterOnsetRootEps*1e6)) microns")
    else
        println("Flutter not found")
    end
    if !isempty(flutterOffsetSpeedsAll)
        iOffLowest = sortperm(flutterOffsetSpeedsAll)[1]
        flutterOffsetSpeed = flutterOffsetSpeedsAll[iOffLowest]
        println("Flutter offset speed = $(round(flutterOffsetSpeed,digits=1)) m/s")
    end

    # Plot V-g-f
    relPath = "/dev/sweptPazy/S0/outputs/PazyWingAirspeedSteps"
    absPath = string(pwd(),relPath)
    mkpath(absPath)
    modeColors = cgrad(:rainbow, nModes, categorical=true)
    ts = 10
    fs = 16
    lw = 2
    gr()
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URangeFlutter), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    for mode in 1:nModes
        plot!(URangeFlutter, modeFrequencies[mode]/(2π), c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URangeFlutter), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    plot!(URangeFlutter, zeros(length(URangeFlutter)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URangeFlutter, modeDampingRatios[mode], c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(plt_Vgf,string(absPath,"/PazyWingAirspeedSteps_",idString,"_Vgf.pdf"))
    gui(plt_Vgf)
end

# --- Dynamic problems setup ---

# Airspeed range for dynamic problems
URangeDyn = collect(floor(flutterOnsetSpeed)+1:0.125:ceil(flutterOffsetSpeed+1))

# Airspeed profile variables
N = length(URangeDyn)
τ₀Range = 1*ones(N-1)  # Steady lower value time
τ₁Range = 1*ones(N-1)  # Linear transition time
τ₂Range = 2*ones(N-1)  # Steady upper value time
τ₂Range[2] = 10

# Time variables
Δt = 2.5e-4
trackingFrequency = round(Int,1e-3/Δt)

# System solver options
maxIter = 50
relTol = 1e-9
NR = create_NewtonRaphson(maximumIterations=maxIter,relativeTolerance=relTol,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# --- Solve initial steady problem ---
println("Solving initial steady problem")

# Create and solve steady problem for initial airspeed
steadyProblem = create_SteadyProblem(model=first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,θ=θ,airspeed=URangeDyn[1],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,tipMass=tipMass,ηtipMass=[0;tipMassPos;0],tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent)))
solve!(steadyProblem)

# --- Solve dynamic problems ---
println("Solving dynamic problems")

# Initialize outputs
dynamicProblem = Array{DynamicProblem}(undef,N-1)
Uvec = Array{Vector{Float64}}(undef,N-1)
t = Array{Vector{Float64}}(undef,N-1)
tipOOP = Array{Vector{Float64}}(undef,N-1)
tipAoA = Array{Vector{Float64}}(undef,N-1)
rootEpsLE = Array{Vector{Float64}}(undef,N-1)
rootEpsTE = Array{Vector{Float64}}(undef,N-1)
timespan = Array{Float64}(undef,N-1)
initialTime = Array{Float64}(undef,N-1)
finalTime = Array{Float64}(undef,N-1)

# Step over airspeeds
for (j,U) in enumerate(URangeDyn[1:end-1])
    # Display progress
    println("Solving for U = $U -> $(URangeDyn[j+1]) m/s")
    # Set initial solution
    x0 = j == 1 ? steadyProblem.x : dynamicProblem[j-1].x
    # Span, initial and final times of simulation
    timespan[j] = τ₀Range[j] + τ₁Range[j] + τ₂Range[j]
    initialTime[j] = j == 1 ? 0.0 : finalTime[j-1]
    finalTime[j] = initialTime[j] + timespan[j]
    # Set airspeed profile
    airspeed = t -> ifelse(t<=initialTime[j]+τ₀Range[j], U, ifelse(t>=initialTime[j]+τ₀Range[j]+τ₁Range[j], URangeDyn[j+1], U + (URangeDyn[j+1]-U)*(t-τ₀Range[j]-initialTime[j])/τ₁Range[j]))
    # Model
    dynModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,θ=θ,airspeed=airspeed,tipMass=tipMass,ηtipMass=[0;tipMassPos;0],tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent))
    # Reference dynamic problem
    refDynamicProblem = j==1 ? nothing : dynamicProblem[j-1]
    # Create and solve dynamic problem
    dynamicProblem[j] = create_DynamicProblem(model=dynModel,initialTime=initialTime[j],finalTime=finalTime[j],Δt=Δt,systemSolver=NR,x0=x0,refDynamicProblem=refDynamicProblem,trackingFrequency=trackingFrequency,saveInitialSolution=false)
    solve!(dynamicProblem[j])
    # Unpack numerical solution
    t[j] = dynamicProblem[j].savedTimeVector
    Uvec[j] = [v[2] for v in dynamicProblem[j].model.v_A.(t[j])]
    tipAoA[j] = [dynamicProblem[j].aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in eachindex(t[j])]
    tipOOP[j] = -[dynamicProblem[j].nodalStatesOverTime[i][nElem].u_n2[1] for i in eachindex(t[j])]
    ϵ11RootDyn = [dynamicProblem[j].compElementalStatesOverTime[i][1].γ[1] for i in eachindex(t[j])]
    κ2RootDyn = [dynamicProblem[j].compElementalStatesOverTime[i][1].κ[2] for i in eachindex(t[j])]
    κ3RootDyn = [dynamicProblem[j].compElementalStatesOverTime[i][1].κ[3] for i in eachindex(t[j])]
    rootEpsLE[j] = (ϵ11RootDyn .- κ2RootDyn*zSG .- κ3RootDyn*ySG_LE)
    rootEpsTE[j] = (ϵ11RootDyn .- κ2RootDyn*zSG .- κ3RootDyn*ySG_TE)
end

# Concatenate dynamic problems outputs
tAll = vcat(t...)
UAll = vcat(Uvec...)
tipAoAAll = vcat(tipAoA...)
tipOOPAll = vcat(tipOOP...)
rootEpsLEAll = vcat(rootEpsLE...)
rootEpsTEAll = vcat(rootEpsTE...)
rootEpsAll = (rootEpsLEAll.+rootEpsTEAll)/2

# Set case ID for dynamic plots
idString = URangeDyn[end] > URangeDyn[1] ? string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_upsteps") : string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name, "_downsteps")

# Set paths
relPath = "/dev/sweptPazy/S0/outputs/PazyWingAirspeedSteps"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
showAirspeedInPlots = true
interactivePlots = false
colors = cgrad(:rainbow, 2, categorical=true)
epsColors = cgrad(:rainbow, 3, categorical=true)
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
savefig(string(absPath,"/PazyWingAirspeedSteps_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=extrema(tAll), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(tAll, tipOOPAll/L*100, lw=lw, c=colors[1], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([tAll[jF]], [tipOOPAll[jF]]/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([tAll[jF], tAll[jF]], [tipOOPAll[1], tipOOPAll[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, tAll[jF]], [tipOOPAll[jF], tipOOPAll[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_tipOOP_twin = twinx(plt_tipOOP)
    plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(tAll), ylims=extrema(UAll))
    plot!(plt_tipOOP_twin, tAll, UAll, lw=lw, c=:black, label=false)
end
display(plt_tipOOP)
savefig(string(absPath,"/PazyWingAirspeedSteps_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=extrema(tAll), tickfont=font(ts), guidefont=font(fs))
plot!(tAll, tipAoAAll*180/π, lw=lw, c=colors[1], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([tAll[jF]], [tipAoAAll[jF]]*180/π, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([tAll[jF], tAll[jF]], [tipAoAAll[1], tipAoAAll[jF]]*180/π, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, tAll[jF]], [tipAoAAll[jF], tipAoAAll[jF]]*180/π, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_tipAOA_twin = twinx(plt_tipAOA)
    plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(tAll), ylims=extrema(UAll))
    plot!(plt_tipAOA_twin, tAll, UAll, lw=lw, c=:black, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/PazyWingAirspeedSteps_",idString,"_tipAoA.pdf"))

# Root axial strains
plt_rootEps = plot(xlabel="Time [s]", ylabel="Root axial strains (\$\\mu\$)", xlims=extrema(tAll), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(tAll, rootEpsLEAll*1e6, lw=lw, c=epsColors[1], alpha=α, label="LE")
plot!(tAll, rootEpsTEAll*1e6, lw=lw, c=epsColors[2], alpha=α, label="TE")
if showAirspeedInPlots && !interactivePlots
    plt_rootEps_twin = twinx(plt_rootEps)
    plot!(plt_rootEps_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(tAll), ylims=extrema(UAll))
    plot!(plt_rootEps_twin, tAll, UAll, lw=lw, c=:black, label=false)
end
display(plt_rootEps)
savefig(string(absPath,"/PazyWingAirspeedSteps_",idString,"_rootEps.pdf"))

# Root strains spectrogram
tBegin = 0
jBegin = findfirst(x-> x>=tBegin, tAll)
tCropped = tAll[jBegin:end] .- tAll[jBegin]
rootEpsCropped = rootEpsAll[jBegin:end]
f_eps, t_eps, S_eps = AeroBeams.spectrogram(tCropped, rootEpsCropped*1e6, window=2^12)
U_eps = UAll[searchsortedfirst.(Ref(tAll), t_eps)]
plt_eps_spectrum = heatmap(U_eps, f_eps, S_eps; ylims=[0,100], title="Root strains spectrogram", xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
plot!([flutterOnsetSpeed, flutterOnsetSpeed], [0, 100], ls=:dash, c=:black, lw=lw, label=false)
plot!([UAll[1], U_eps[end]], [flutterOnsetFreq, flutterOnsetFreq]/(2π), ls=:dash, c=:black, lw=lw, label=false)
plot!([UAll[1], U_eps[end]], 2*[flutterOnsetFreq, flutterOnsetFreq]/(2π), ls=:dash, c=:black, lw=lw, label=false)
plot!([UAll[1], U_eps[end]], 3*[flutterOnsetFreq, flutterOnsetFreq]/(2π), ls=:dash, c=:black, lw=lw, label=false)
display(plt_eps_spectrum)
savefig(string(absPath,"/PazyWingAirspeedSteps_",idString,"_rootEpsSpectrogram.pdf"))

println("Finished PazyWingAirspeedSteps.jl")