using AeroBeams, LinearInterpolations, Plots, ColorSchemes, JLD2

# Select tip mass configuration (set as XX_mY_oZ, where XX is either LE or TE, Y is the mass in grams and Z is the offset in cm)
tipMassConfig = "LE_m0_o0"

# Airspeed sweep type (choose between "upsweep" and "downsweep")
UsweepType = "upsweep"
@assert UsweepType in ["upsweep","downsweep"]

# Root pitch angle range
θRange = π/180*vcat(8.5)

# Airspeed step in dynamic problems
ΔU = 0.5

# Fixed time variables
Δt = 2.5e-4
tSim = 25
trackingFrequency = round(Int,5e-4/Δt)

# Airspeed range for flutter problems
URangeFlutter = collect(20:0.25:90)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "Exponential"

# Aerodynamic solver
aeroSolver = BLi()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Tip impulse configuration to induce LCO
F₀ = 1
ω = 4*2π
τ = 2π/ω
t₀ = 0.1
F = t -> ifelse.(t.<=t₀,0.0,ifelse.(t.<=t₀+τ/2,F₀*sin.(ω*(t.-(t₀+τ))),0.0))

# System solver options for dynamic problems
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(maximumIterations=maxIter,relativeTolerance=relTol,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2,allowAdvanceThroughUnconvergedAeroStates=aeroSolver.name=="BLi")

# Number of modes for flutter computation
nModes = 3

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
ySG_LE = 50e-3/2    # y-position: on LE
ySG_TE = -50e-3/2   # y-position: on TE
zSG = 2.25e-3/2     # z-position: on top

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Set tip mass configuration variables
m = match(r"^(LE|TE)_m(\d+(?:\.\d+)?)_o(\d+(?:\.\d+)?)$", tipMassConfig)
@assert m !== nothing "Invalid tipMassConfig string"
pos = m.captures[1]
tipMass = 1e-3*parse(Float64, m.captures[2])
tipMassPosOffset = 1e-2*parse(Float64, m.captures[3])
tipMassPos = pos == "LE" ? chord*normSparPos + tipMassPosOffset : -(chord*(1-normSparPos) + tipMassPosOffset)

# Initialize flutter outputs
flutterProblem = Array{EigenProblem}(undef,length(θRange),length(URangeFlutter))
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URangeFlutter))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URangeFlutter))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URangeFlutter))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URangeFlutter))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URangeFlutter))
tipOOP_f = Array{Float64}(undef,length(θRange),length(URangeFlutter))
tipAoA_f = Array{Float64}(undef,length(θRange),length(URangeFlutter))
rootEps_f = Array{Float64}(undef,length(θRange),length(URangeFlutter))
flutterOnsetSpeedOfMode = [fill(0.0,0) for _ in eachindex(θRange), _ in 1:nModes]
flutterOffsetSpeedOfMode = [fill(0.0,0) for _ in eachindex(θRange), _ in 1:nModes]
flutterOnsetFreqOfMode = [fill(0.0,0) for _ in eachindex(θRange), _ in 1:nModes]
flutterOnsetTipOOPOfMode = [fill(0.0,0) for _ in eachindex(θRange), _ in 1:nModes]
flutterOnsetTipAoAOfMode = [fill(0.0,0) for _ in eachindex(θRange), _ in 1:nModes]
flutterOnsetRootEpsOfMode = [fill(0.0,0) for _ in eachindex(θRange), _ in 1:nModes]
flutterOnsetSpeed = Array{Float64}(undef,length(θRange))
flutterOnsetFreq = Array{Float64}(undef,length(θRange))
flutterOnsetTipOOP = Array{Float64}(undef,length(θRange))
flutterOnsetTipAoA = Array{Float64}(undef,length(θRange))
flutterOnsetRootEps = Array{Float64}(undef,length(θRange))
flutterOffsetSpeed = Array{Float64}(undef,length(θRange))

# Set paths
relPath = "/dev/sweptPazy/S0/outputs/PazyWingPitchRangeURangeLCO"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
epsColors = cgrad(:rainbow, 2, categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 10
lw = 2
α = 0.5
spectrogramWindow = Int(2^12)
gr()

# Sweep root pitch angle
for (i,θ) in enumerate(θRange)
    # Case ID
    idStringFlutter = string(tipMassConfig, "_Theta", round(θ*180/π,digits=1), "_", tipLossType, "_", aeroSolver.name)
    # --- Solve flutter problem ---
    println("Solving flutter problem for θ = $(round(θ*180/π,digits=1)) deg")
    # Sweep airspeed
    for (j,U) in enumerate(URangeFlutter)
        # Model
        flutterModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,θ=θ,airspeed=U,tipMass=tipMass,ηtipMass=[0;tipMassPos;0]))
        # Set initial guess solution as previously available one
        x0 = j == 1 ? zeros(0) : flutterProblem[i,j-1].x
        # Create and solve problem
        flutterProblem[i,j] = create_EigenProblem(model=flutterModel,nModes=nModes,frequencyFilterLimits=[2π,Inf],x0=x0)
        solve!(flutterProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = flutterProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = flutterProblem[i,j].dampingsOscillatory
        untrackedEigenvectors[i,j] = flutterProblem[i,j].eigenvectorsOscillatoryCplx
        # Kinematics
        tipOOP_f[i,j] = -flutterProblem[i,j].nodalStatesOverσ[end][nElem].u_n2[1]
        tipAoA_f[i,j] = flutterProblem[i,j].aeroVariablesOverσ[end][nElem].flowAnglesAndRates.αₑ
        ϵ11Root_f = flutterProblem[i,j].compElementalStatesOverσ[end][1].γ[1]
        κ2Root_f = flutterProblem[i,j].compElementalStatesOverσ[end][1].κ[2]
        rootEps_f[i,j] = (ϵ11Root_f .- κ2Root_f*zSG)
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(URangeFlutter,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and damping ratios by mode
    modeFrequencies = Array{Vector{Float64}}(undef,length(θRange),nModes)
    modeDampings = Array{Vector{Float64}}(undef,length(θRange),nModes)
    modeDampingRatios = Array{Vector{Float64}}(undef,length(θRange),nModes)
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URangeFlutter)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URangeFlutter)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
    # Compute flutter variables
    for mode in 1:nModes
        jOnset = findall(j -> modeDampings[i,mode][j] < 0 && modeDampings[i,mode][j+1] > 0, 1:length(URangeFlutter)-1)
        jOffset = findall(j -> modeDampings[i,mode][j] > 0 && modeDampings[i,mode][j+1] < 0, 1:length(URangeFlutter)-1)
        if !isempty(jOnset)
            for jO in jOnset
                push!(flutterOnsetSpeedOfMode[i,mode],interpolate(modeDampings[i,mode][jO:jO+1],URangeFlutter[jO:jO+1], 0))
                push!(flutterOnsetFreqOfMode[i,mode],interpolate(modeDampings[i,mode][jO:jO+1],modeFrequencies[i,mode][jO:jO+1], 0))
                push!(flutterOnsetTipOOPOfMode[i,mode], tipOOP_f[i,jO])
                push!(flutterOnsetTipAoAOfMode[i,mode], tipAoA_f[i,jO])
                push!(flutterOnsetRootEpsOfMode[i,mode], rootEps_f[i,jO])
            end
        end
        if !(isempty(jOffset) || isempty(jOnset))
            for jO in jOffset
                push!(flutterOffsetSpeedOfMode[i,mode],interpolate(-modeDampings[i,mode][jO:jO+1],URangeFlutter[jO:jO+1],0))
            end
        end
    end
    flutterOnsetSpeedsAll = vcat(filter(!isempty,flutterOnsetSpeedOfMode[i,:])...)
    flutterOffsetSpeedsAll = vcat(filter(!isempty,flutterOffsetSpeedOfMode[i,:])...)
    flutterOnsetFreqAll = vcat(filter(!isempty,flutterOnsetFreqOfMode[i,:])...)
    flutterOnsetTipOOPAll = vcat(filter(!isempty,flutterOnsetTipOOPOfMode[i,:])...)
    flutterOnsetTipAoAAll = vcat(filter(!isempty,flutterOnsetTipAoAOfMode[i,:])...)
    flutterOnsetRootEpsAll = vcat(filter(!isempty,flutterOnsetRootEpsOfMode[i,:])...)
    if !isempty(flutterOnsetSpeedsAll)
        iOnLowest = sortperm(flutterOnsetSpeedsAll)[1]
        flutterOnsetSpeed[i] = flutterOnsetSpeedsAll[iOnLowest]
        flutterOnsetFreq[i] = flutterOnsetFreqAll[iOnLowest]
        flutterOnsetTipOOP[i] = flutterOnsetTipOOPAll[iOnLowest]
        flutterOnsetTipAoA[i] = flutterOnsetTipAoAAll[iOnLowest]
        flutterOnsetRootEps[i] = flutterOnsetRootEpsAll[iOnLowest]
        println("Flutter onset speed = $(round(flutterOnsetSpeed[i],digits=1)) m/s, frequency = $(round(flutterOnsetFreq[i]/(2*π),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP[i]/L*100,digits=1))% semispan, tip AoA = $(round(flutterOnsetTipAoA[i]*180/π,digits=1)) deg, root strains = $(round(Int,flutterOnsetRootEps[i]*1e6)) microns")
    else
        println("Flutter not found")
    end
    if !isempty(flutterOffsetSpeedsAll)
        iOffLowest = sortperm(flutterOffsetSpeedsAll)[1]
        flutterOffsetSpeed[i] = flutterOffsetSpeedsAll[iOffLowest]
        println("Flutter offset speed = $(round(flutterOffsetSpeed[i],digits=1)) m/s")
    end
    # Plot V-g-f
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URangeFlutter), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    for mode in 1:nModes
        plot!(URangeFlutter, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URangeFlutter), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    plot!(URangeFlutter, zeros(length(URangeFlutter)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URangeFlutter, modeDampingRatios[i,mode], c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(plt_Vgf,string(absPath,"/PazyWingAirspeedSteps_",idStringFlutter,"_Vgf.pdf"))
    gui(plt_Vgf)

    # --- Dynamic problems setup ---
    # Airspeed range for dynamic problems
    URangeDyn = UsweepType == "upsweep" ? collect(39.5:ΔU:40.5) : collect(ceil(flutterOffsetSpeed[i])+ΔU:-ΔU:floor(flutterOnsetSpeed[i])-ΔU)
    # Set tip impulse (on a dummy beam, updated later on model creation)
    impulse = create_BC(name="impulse",beam=create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1)]),node=nElem+1,types=["F1A"],values=[t->F(t)])
    # Step over airspeeds
    for (j,U) in enumerate(URangeDyn)
        # Display progress
        println("Solving $UsweepType dynamic problem for θ = $(round(θ*180/π,digits=1)) deg, U = $U m/s")
        # Case ID
        idStringDyn = string(idStringFlutter, "_", UsweepType, "_U", round(U,digits=3))
        # Model
        dynModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,θ=θ,airspeed=U,tipMass=tipMass,ηtipMass=[0;tipMassPos;0],additionalBCs=[impulse]))
        # Create and solve steady problem for initial solution
        steadyProblem = create_SteadyProblem(model=dynModel)
        solve!(steadyProblem)
        # Create and solve dynamic problem
        dynamicProblem = create_DynamicProblem(model=dynModel,finalTime=tSim,Δt=Δt,systemSolver=NR,x0=steadyProblem.x,trackingFrequency=trackingFrequency)
        solve!(dynamicProblem)
        # Unpack numerical solution
        t = dynamicProblem.savedTimeVector
        tip_p = [dynamicProblem.nodalStatesOverTime[k][nElem].p_n2_b for k in eachindex(t)]
        tipTwist = [asin((first(rotation_tensor_WM(tip_p[k]))*AeroBeams.a2)[3]) for k in eachindex(t)]
        tipAoA = [dynamicProblem.aeroVariablesOverTime[k][nElem].flowAnglesAndRates.αₑ for k in eachindex(t)]
        tipOOP = -[dynamicProblem.nodalStatesOverTime[k][nElem].u_n2[1] for k in eachindex(t)]
        ϵ11RootDyn = [dynamicProblem.compElementalStatesOverTime[k][1].γ[1] for k in eachindex(t)]
        κ2RootDyn = [dynamicProblem.compElementalStatesOverTime[k][1].κ[2] for k in eachindex(t)]
        κ3RootDyn = [dynamicProblem.compElementalStatesOverTime[k][1].κ[3] for k in eachindex(t)]
        rootEpsLE = (ϵ11RootDyn .- κ2RootDyn*zSG .- κ3RootDyn*ySG_LE)
        rootEpsTE = (ϵ11RootDyn .- κ2RootDyn*zSG .- κ3RootDyn*ySG_TE)
        # Save data
        @save absPath*"/"*idStringDyn*"_t.jld2" t
        @save absPath*"/"*idStringDyn*"_tipOOP.jld2" tipOOP
        @save absPath*"/"*idStringDyn*"_tipAoA.jld2" tipAoA
        @save absPath*"/"*idStringDyn*"_tipTwist.jld2" tipTwist
        @save absPath*"/"*idStringDyn*"_rootEpsLE.jld2" rootEpsLE
        @save absPath*"/"*idStringDyn*"_rootEpsTE.jld2" rootEpsTE
        # Tip OOP displacement
        title = "θ = $(round(θ*180/π,digits=1)) deg, U = $U m/s, $UsweepType"
        plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
        plot!(t, tipOOP/L*100, lw=lw, c=:black, label=false)
        display(plt_tipOOP)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_tipOOP.pdf"))
        # Tip AoA
        plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs))
        plot!(t, tipAoA*180/π, lw=lw, c=:black, label=false)
        display(plt_tipAOA)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_tipAoA.pdf"))
        # Tip twist
        plt_tipTwist = plot(xlabel="Time [s]", ylabel="Tip twist [deg]", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs))
        plot!(t, tipTwist*180/π, lw=lw, c=:black, label=false)
        display(plt_tipTwist)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_tipTwist.pdf"))
        # Root axial strains
        plt_rootEps = plot(xlabel="Time [s]", ylabel="Root axial strains (\$\\mu\$)", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
        plot!(t, rootEpsLE*1e6, lw=lw, c=epsColors[1], alpha=α, label="LE")
        plot!(t, rootEpsTE*1e6, lw=lw, c=epsColors[2], alpha=α, label="TE")
        display(plt_rootEps)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_rootEps.pdf"))
        # Root strains spectrogram
        rootEps = (rootEpsLE.+rootEpsTE)/2
        f_eps, t_eps, S_eps = spectrogram(t, rootEps*1e6, window=spectrogramWindow)
        plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title=title, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
        display(plt_eps_spectrum)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_rootEpsSpectrogram.pdf"))
        GC.gc()
    end
end

println("Finished PazyWingPitchRangeURangeLCO.jl")