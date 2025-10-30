using AeroBeams, LinearInterpolations, Plots, ColorSchemes

# Select tip mass configuration (set as XX_mY_oZ, where XX is either LE or TE, Y is the mass in grams and Z is the offset in cm)
tipMassConfig = "LE_m0_o0"

# Flag to solve flutter problem
solveFlutter = true

# Root pitch angle offset (assumed)
θoffset = 0.0*π/180

# Root pitch angle
θ = 7π/180 + θoffset

# Profile for airspeed sweep
τ = 1
Umin = 42
Umax = 42.5
U = t -> ifelse(t<=τ, Umin + (Umax-Umin)*t/τ, Umax)

# Airspeed range for flutter problem
URange = collect(20:0.25:90)

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

# Time variables
Δt = 5e-4
tf = 20
trackingFrequency = round(Int,1e-3/Δt)

# Number of modes for flutter computation
nModes = 3

# Set system solver options for dynamic problem
maxIter = 50
relTol = 1e-9
NR = create_NewtonRaphson(maximumIterations=maxIter,relativeTolerance=relTol,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
ySG_LE = 50e-3/2    # y-position: on LE
ySG_TE = -50e-3/2   # y-position: on TE
ySG = 0             # average over LE and TE
zSG = 2.25e-3/2     # z-position: on top

# Case ID
idString = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_Umax", Umax, "_", tipLossType, "_", aeroSolver.name)

# --- Solve flutter problem ---
if solveFlutter
    println("Solving flutter problem")

    # Initialize outputs
    flutterProblem = Array{EigenProblem}(undef,length(URange))
    untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
    untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
    untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
    tipOOP_f = Array{Float64}(undef,length(URange))
    tipAoA_f = Array{Float64}(undef,length(URange))
    rootEps_f = Array{Float64}(undef,length(URange))

    # Sweep airspeed
    for (j,U) in enumerate(URange)
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
        κ3Root_f = flutterProblem[j].compElementalStatesOverσ[end][1].κ[3]
        rootEps_f[j] = (ϵ11Root_f .- κ2Root_f*zSG .- κ3Root_f*ySG)
    end

    # Apply mode tracking
    freqs,damps,_ = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

    # Separate frequencies and damping ratios by mode
    modeFrequencies = Array{Vector{Float64}}(undef,nModes)
    modeDampings = Array{Vector{Float64}}(undef,nModes)
    modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
    for mode in 1:nModes
        modeFrequencies[mode] = [freqs[j][mode] for j in eachindex(URange)]
        modeDampings[mode] = [damps[j][mode] for j in eachindex(URange)]
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
        jOnset = findall(j -> modeDampings[mode][j] < 0 && modeDampings[mode][j+1] > 0, 1:length(URange)-1)
        jOffset = findall(j -> modeDampings[mode][j] > 0 && modeDampings[mode][j+1] < 0, 1:length(URange)-1)
        if !isempty(jOnset)
            for jO in jOnset
                push!(flutterOnsetSpeedOfMode[mode],interpolate(modeDampings[mode][jO:jO+1],URange[jO:jO+1], 0))
                push!(flutterOnsetFreqOfMode[mode],interpolate(modeDampings[mode][jO:jO+1],modeFrequencies[mode][jO:jO+1], 0))
                push!(flutterOnsetTipOOPOfMode[mode], tipOOP_f[jO])
                push!(flutterOnsetTipAoAOfMode[mode], tipAoA_f[jO])
                push!(flutterOnsetRootEpsOfMode[mode], rootEps_f[jO])
            end
        end
        if !(isempty(jOffset) || isempty(jOnset))
            for jO in jOffset
                push!(flutterOffsetSpeedOfMode[mode],interpolate(-modeDampings[mode][jO:jO+1],URange[jO:jO+1],0))
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
    relPath = "/dev/sweptPazy/S0/outputs/PazyWingAirspeedSweep"
    absPath = string(pwd(),relPath)
    mkpath(absPath)
    modeColors = cgrad(:rainbow, nModes, categorical=true)
    ts = 10
    fs = 16
    lw = 2
    gr()
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    for mode in 1:nModes
        plot!(URange, modeFrequencies[mode]/(2π), c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(plt_Vgf,string(absPath,"/PazyWingAirspeedSweep_",idString,"_Vgf.pdf"))
    gui(plt_Vgf)
end

# --- Solve dynamic problem ---
println("Solving dynamic problem")

# Model
dynModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,θ=θ,airspeed=U,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,tipMass=tipMass,ηtipMass=[0;tipMassPos;0],tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent))

# Create and solve steady problem for initial airspeed
steadyProblem = create_SteadyProblem(model=dynModel)
solve!(steadyProblem)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=dynModel,finalTime=tf,Δt=Δt,systemSolver=NR,x0=steadyProblem.x,trackingFrequency=trackingFrequency)
solve!(dynamicProblem)

# Unpack numerical solution
t = dynamicProblem.savedTimeVector
tipAoA = [dynamicProblem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in eachindex(t)]
tipOOP = -[dynamicProblem.nodalStatesOverTime[i][nElem].u_n2[1] for i in eachindex(t)]
ϵ11Root = [dynamicProblem.compElementalStatesOverTime[i][1].γ[1] for i in eachindex(t)]
κ2Root = [dynamicProblem.compElementalStatesOverTime[i][1].κ[2] for i in eachindex(t)]
κ3Root = [dynamicProblem.compElementalStatesOverTime[i][1].κ[3] for i in eachindex(t)]
rootEps = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG)
rootEpsLE = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG_LE)
rootEpsTE = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG_TE)

# Set paths
relPath = "/dev/sweptPazy/S0/outputs/PazyWingAirspeedSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(dynamicProblem,refBasis="A",plotFrequency=round(Int,5e-2/Δt),plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingAirspeedSweep_",idString,"_deformation.gif"),displayProgress=true)

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
savefig(string(absPath,"/PazyWingAirspeedSweep_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, tipOOP/L*100, lw=lw, c=colors[1], label="AeroBeams")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=colors[2], label="Revivo & Raveh (2025)")
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], [tipOOP[jF]]/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], [tipOOP[1], tipOOP[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, t[jF]], [tipOOP[jF], tipOOP[jF]]/L*100, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_tipOOP_twin = twinx(plt_tipOOP)
    plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_tipOOP_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipOOP)
savefig(string(absPath,"/PazyWingAirspeedSweep_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, tipAoA*180/π, lw=lw, c=colors[1], label=false)
plot!([NaN], [NaN], lw=lw, ls=:dash, c=colors[2], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], [tipAoA[jF]]*180/π, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], [tipAoA[1], tipAoA[jF]]*180/π, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, t[jF]], [tipAoA[jF], tipAoA[jF]]*180/π, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_tipAOA_twin = twinx(plt_tipAOA)
    plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_tipAOA_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/PazyWingAirspeedSweep_",idString,"_tipAoA.pdf"))

# Root axial strains
plt_rootEps = plot(xlabel="Time [s]", ylabel="Root axial strains (\$\\mu\$)", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, rootEps*1e6, lw=lw, c=epsColors[1], alpha=0.7, label="Avg.")
plot!(t, rootEpsLE*1e6, lw=lw, c=epsColors[2], alpha=0.7, label="LE")
plot!(t, rootEpsTE*1e6, lw=lw, c=epsColors[3], alpha=0.7, label="TE")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=epsColors[2], label=false)
if solveFlutter && !isnothing(jF)
    scatter!([t[jF]], [rootEps[jF]]*1e6, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
    plot!([t[jF], t[jF]], [rootEps[1], rootEps[jF]]*1e6, ls=:dash, c=:black, lw=lw, label=false)
    plot!([0, t[jF]], [rootEps[jF], rootEps[jF]]*1e6, ls=:dash, c=:black, lw=lw, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_rootEps_twin = twinx(plt_rootEps)
    plot!(plt_rootEps_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_rootEps_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_rootEps)
savefig(string(absPath,"/PazyWingAirspeedSweep_",idString,"_rootEps.pdf"))

# Root strains spectrogram
tBegin = 0
jBegin = findfirst(x-> x>=tBegin, t)
tCropped = t[jBegin:end] .- t[jBegin]
rootEpsCropped = rootEps[jBegin:end]
f_eps, t_eps, S_eps = AeroBeams.spectrogram(tCropped, rootEpsCropped*1e6, window=2^12)
U_vec = U.(tCropped)
U_eps = U_vec[searchsortedfirst.(Ref(tCropped), t_eps)]
plt_eps_spectrum = heatmap(U_eps, f_eps, S_eps; ylims=[0,100], title="Root strains spectrogram", xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
plot!([U_vec[1], U_vec[end]], [flutterOnsetFreq, flutterOnsetFreq]/(2π), ls=:dash, c=:black, lw=lw, label=false)
plot!([U_vec[1], U_vec[end]], 2*[flutterOnsetFreq, flutterOnsetFreq]/(2π), ls=:dash, c=:black, lw=lw, label=false)
plot!([U_vec[1], U_vec[end]], 3*[flutterOnsetFreq, flutterOnsetFreq]/(2π), ls=:dash, c=:black, lw=lw, label=false)
display(plt_eps_spectrum)
savefig(string(absPath,"/PazyWingAirspeedSweep_",idString,"_rootEpsSpectrogram.pdf"))

println("Finished PazyWingAirspeedSweep.jl")