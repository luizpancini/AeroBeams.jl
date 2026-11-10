using AeroBeams, JLD2

# Select tip mass configuration (set as XX_mY_oZ, where XX is either LE or TE, Y is the mass in grams and Z is the offset in cm)
tipMassConfig = "LE_m0_o0"

# Airspeed sweep type (choose between "upsweep" and "downsweep")
UsweepType = "upsweep"
@assert UsweepType in ["upsweep","downsweep"]

# Root pitch angle range [rad]
θRange = deg2rad.(vcat(7))

# Airspeed step in dynamic problems [m/s]
ΔU = 0.5

# Fixed time variables
Δt = 2.5e-4
tSim = 25
trackingFrequency = round(Int,5e-4/Δt)

# Airspeed range for flutter problems [m/s]
URangeFlutter = collect(20:0.25:90)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "Exponential"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Tip impulse configuration to induce LCO
F₀ = 1
ω = Hz2rads(4)
τ = 2π/ω
t₀ = 0.1
F = t -> ifelse.(t.<=t₀,0.0,ifelse.(t.<=t₀+τ/2,F₀*sin.(ω*(t.-(t₀+τ))),0.0))

# System solver options for dynamic problems
systemSolver = create_NLsolver()

# Number of modes for flutter computation
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Strain gauge coordinates (the spar cs is 60 x 2.25 mm)
SGelem = 3          # element number where the strain gauge is located (see Fig. 12 of https://doi.org/10.2514/1.J060621: FBG^1_3 is @ 127 mm spanwise, but wing root is at 32 mm -> actual spanwise position is 127 - 32 = 95 mm which is at the midpoint of element 3)
x2SG = 50e-3/2      # chordwise position: close to LE
x3SG = 2.25e-3/2    # height: on top

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
eigensols = Array{EigenSolution}(undef,length(θRange),length(URangeFlutter))
tipOOP_f = Array{Float64}(undef,length(θRange),length(URangeFlutter))
tipAoA_f = Array{Float64}(undef,length(θRange),length(URangeFlutter))
epsSG_f = Array{Float64}(undef,length(θRange),length(URangeFlutter))
flutterOnsetSpeed = Array{Float64}(undef,length(θRange))
flutterOnsetFreq = Array{Float64}(undef,length(θRange))
flutterOnsetTipOOP = Array{Float64}(undef,length(θRange))
flutterOnsetTipAoA = Array{Float64}(undef,length(θRange))
flutterOnsetEps = Array{Float64}(undef,length(θRange))
flutterOffsetSpeed = Array{Float64}(undef,length(θRange))

## Plots
# Set paths
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingPitchRangeURangeLCO"
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
    idStringFlutter = string(tipMassConfig, "_Theta", round(rad2deg(θ),digits=1), "_", tipLossType, "_", aeroSolver.name)
    # --- Solve flutter problem ---
    println("Solving flutter problem for θ = $(round(rad2deg(θ),digits=1)) deg")
    # Sweep airspeed
    for (j,U) in enumerate(URangeFlutter)
        # Model
        flutterModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,pitchAngle=θ,airspeed=U,tipMasses=[tipMass],ηtipMasses=[[0;tipMassPos;0]],tipMassesInertias=[zeros(3,3)]))
        # Set initial guess solution as previously available one
        x0 = j == 1 ? zeros(0) : flutterProblem[i,j-1].systemArrays.x
        # Create and solve problem
        flutterProblem[i,j] = create_EigenProblem(model=flutterModel,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(flutterProblem[i,j])
        # Eigensolution
        eigensols[i,j] = flutterProblem[i,j].eigensolution
        # Kinematics
        states = flutterProblem[i,j].loadstepHistory.states
        tipOOP_f[i,j] = -states[end][nElem].nodal.u2[1]
        tipAoA_f[i,j] = flutterProblem[i,j].loadstepHistory.aeroVariables[end][nElem].unsteadyKin.αₑ
        ϵ11SG_f = states[end][SGelem].midpoint.γ[1]
        κ2SG_f = states[end][SGelem].midpoint.κ[2]
        κ3SG_f = states[end][SGelem].midpoint.κ[3]
        epsSG_f[i,j] = (ϵ11SG_f .- κ2SG_f*x3SG .- κ3SG_f*x2SG)
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=2π, maxTracks=nModes)
    # Damping ratios
    dampingRatios = dampsOsc./freqsOsc
    # Flutter onset and offset speeds, frequencies and kinematics
    flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterOnsetTipOOPs,flutterOnsetTipAoAs,flutterOnsetEpss,flutterOffsetSpeeds,_ = find_flutter(URangeFlutter,freqsOsc,dampsOsc,tipOOP_f[i,:],tipAoA_f[i,:],epsSG_f[i,:])
    if !isempty(flutterOnsetSpeeds)
        flutterOnsetSpeed[i] = flutterOnsetSpeeds[1]
        flutterOnsetFreq[i] = flutterOnsetFreqs[1]
        flutterOnsetTipOOP[i] = flutterOnsetTipOOPs[1]
        flutterOnsetTipAoA[i] = flutterOnsetTipAoAs[1]
        flutterOnsetEps[i] = flutterOnsetEpss[1]
        println("Flutter onset speed = $(round(flutterOnsetSpeed[i],digits=1)) m/s, frequency = $(round(rads2Hz(flutterOnsetFreq[i]),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP[i]/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOnsetTipAoA[i]),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOnsetEps[i]*1e6)) microns")
    else
        println("Flutter not found")
    end
    if !isempty(flutterOffsetSpeeds)
        flutterOffsetSpeed[i] = flutterOffsetSpeeds[1]
        println("Flutter offset speed = $(round(flutterOffsetSpeed[i],digits=1)) m/s")
    end
    # Plot V-g-f
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
    savefig(plt_Vgf,string(absPath,"/PazyWingAirspeedSteps_",idStringFlutter,"_Vgf.pdf"))

    # --- Dynamic problems setup ---
    # Airspeed range for dynamic problems [m/s]
    range = collect(floor(flutterOnsetSpeed[i]/ΔU)*ΔU:ΔU:ceil((flutterOffsetSpeed[i]+ΔU)/ΔU)*ΔU)
    URangeDyn = UsweepType == "upsweep" ?
                range :
                reverse(range)
    # Set tip impulse (on a dummy beam, updated later on model creation)
    impulse = create_BC(name="impulse",beam=create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1)]),node=nElem+1,conditions=Dict(:F1A=>t->F(t)))
    # Step over airspeeds
    for (j,U) in enumerate(URangeDyn)
        # Display progress
        println("Solving $UsweepType dynamic problem for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Case ID
        idStringDyn = string(idStringFlutter, "_", UsweepType, "_U", round(U,digits=3))
        # Model
        dynModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,pitchAngle=θ,airspeed=U,tipMasses=[tipMass],ηtipMasses=[[0;tipMassPos;0]],tipMassesInertias=[zeros(3,3)],additionalBCs=[impulse]))
        # Create and solve steady problem for initial solution
        steadyProblem = create_SteadyProblem(model=dynModel)
        solve!(steadyProblem)
        # Create and solve dynamic problem
        dynamicProblem = create_DynamicProblem(model=dynModel,finalTime=tSim,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x,trackingFrequency=trackingFrequency)
        solve!(dynamicProblem)
        # Unpack numerical solution
        sol = dynamicProblem.solutionHistory
        t = sol.t
        states = sol.states
        aeroVars = sol.aeroVariables
        tipTwist = [twist_angle(states[k][nElem].nodal.p2_b, outputUnit=:rad) for k in eachindex(t)]
        tipAoA = [aeroVars[k][nElem].unsteadyKin.αₑ for k in eachindex(t)]
        tipOOP = -[states[k][nElem].nodal.u2[1] for k in eachindex(t)]
        ϵ11SG = [states[k][SGelem].midpoint.γ[1] for k in eachindex(t)]
        κ2SG = [states[k][SGelem].midpoint.κ[2] for k in eachindex(t)]
        κ3SG = [states[k][SGelem].midpoint.κ[3] for k in eachindex(t)]
        epsSG = (ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG)
        # Save data
        @save absPath*"/"*idStringDyn*"_t.jld2" t
        @save absPath*"/"*idStringDyn*"_tipOOP.jld2" tipOOP
        @save absPath*"/"*idStringDyn*"_tipAoA.jld2" tipAoA
        @save absPath*"/"*idStringDyn*"_tipTwist.jld2" tipTwist
        @save absPath*"/"*idStringDyn*"_epsSG.jld2" epsSG
        # Tip OOP displacement
        title = "θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s, $UsweepType"
        plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
        plot!(t, tipOOP/L*100, lw=lw, c=:black, label=false)
        display(plt_tipOOP)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_tipOOP.pdf"))
        # Tip AoA
        plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs))
        plot!(t, rad2deg.(tipAoA), lw=lw, c=:black, label=false)
        display(plt_tipAOA)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_tipAoA.pdf"))
        # Tip twist
        plt_tipTwist = plot(xlabel="Time [s]", ylabel="Tip twist [deg]", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs))
        plot!(t, rad2deg.(tipTwist), lw=lw, c=:black, label=false)
        display(plt_tipTwist)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_tipTwist.pdf"))
        # 1/6-span strains
        plt_eps = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\\epsilon\$)", title=title, xlims=[0,tSim], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
        plot!(t, epsSG*1e6, lw=lw, c=epsColors[1], alpha=α, label=false)
        display(plt_eps)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_eps.pdf"))
        # 1/6-span strains spectrogram
        f_eps, t_eps, S_eps = spectrogram(t, epsSG*1e6, window=spectrogramWindow)
        plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title=title, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
        display(plt_eps_spectrum)
        savefig(string(absPath,"/PazyWingPitchRangeURangeLCO_",idStringDyn,"_epsSpectrogram.pdf"))
        GC.gc()
    end
end

println("Finished PazyWingPitchRangeURangeLCO.jl")
