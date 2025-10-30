using AeroBeams, LinearInterpolations, Plots, ColorSchemes

# Tip mass configurations (assumed tipMassPosOffset ahead/behind the LE/TE)
tipMassConfigs = ["LE", "TE"]
tipMasses = [11e-3, 11e-3]
tipMassPosOffsets = [0.05, 0.05]

# Root pitch angle range
θRange = π/180*vcat(1,3,5,7,8.5,10)

# Airspeed range
URange = collect(30:1:85)

# Sweep angle [rad]
Λ = 10*π/180

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Number of modes
nModes = 3

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),nModes)
tipOOP = Array{Float64}(undef,length(tipMassConfigs),length(θRange),length(URange))

# Sweep tip mass configurations
for (c,config) in enumerate(tipMassConfigs)
    # Set current tip mass position along the chord
    tipMassPos = config == "LE" ? (chord*normSparPos + tipMassPosOffsets[c]) : (- (chord*(1-normSparPos) + tipMassPosOffsets[c]))
    # Sweep root pitch angle
    for (i,θ) in enumerate(θRange)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for tip mass at $config, θ = $(round(θ*180/pi,digits=1)) deg, U = $U m/s")
            # Model
            model = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMasses[c],ηtipMass=[0;tipMassPos;0]))
            # Create and solve problem
            problem = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[2π,Inf])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[c,i,j] = problem.frequenciesOscillatory
            untrackedDamps[c,i,j] = problem.dampingsOscillatory
            untrackedEigenvectors[c,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP deflection
            tipOOP[c,i,j] = -problem.nodalStatesOverσ[end][nElem].u_n2[1]
        end
        # Apply mode tracking
        freqs[c,i,:],damps[c,i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[c,i,:],untrackedDamps[c,i,:],untrackedEigenvectors[c,i,:])
        # Separate frequencies and damping ratios by mode
        for mode in 1:nModes
            modeFrequencies[c,i,mode] = [freqs[c,i,j][mode] for j in eachindex(URange)]
            modeDampings[c,i,mode] = [damps[c,i,j][mode] for j in eachindex(URange)]
            modeDampingRatios[c,i,mode] = modeDampings[c,i,mode]./modeFrequencies[c,i,mode]
        end
    end
end

# Compute flutter variables
flutterOnsetMode = [fill(0,0) for _ in 1:length(tipMassConfigs), _ in 1:length(θRange)]
flutterOnsetSpeedOfMode = [fill(0.0,0) for _ in 1:length(tipMassConfigs), _ in 1:length(θRange), _ in 1:nModes]
flutterOnsetFreqOfMode = [fill(0.0,0) for _ in 1:length(tipMassConfigs), _ in 1:length(θRange), _ in 1:nModes]
flutterOnsetTipOOPOfMode = [fill(0.0,0) for _ in 1:length(tipMassConfigs), _ in 1:length(θRange), _ in 1:nModes]
flutterOnsetSpeed = [Inf for _ in 1:length(tipMassConfigs), _ in 1:length(θRange)]
flutterOnsetFreq = [NaN for _ in 1:length(tipMassConfigs), _ in 1:length(θRange)]
flutterOnsetTipOOP = [NaN for _ in 1:length(tipMassConfigs), _ in 1:length(θRange)]
flutterOnsetSpeedsAll = [fill(0.0,0) for _ in 1:length(tipMassConfigs), _ in 1:length(θRange)]
flutterOnsetFreqsAll = [fill(0.0,0) for _ in 1:length(tipMassConfigs), _ in 1:length(θRange)]
for (c,config) in enumerate(tipMassConfigs)
    for (i,θ) in enumerate(θRange)
        # Flutter onset/offset data of each mode
        for mode in 1:nModes
            jOnset = findall(j -> modeDampings[c,i,mode][j] < 0 && modeDampings[c,i,mode][j+1] > 0, 1:length(URange)-1)
            if !isempty(jOnset)
                for jO in jOnset
                    push!(flutterOnsetMode[c,i],mode)
                    push!(flutterOnsetSpeedOfMode[c,i,mode],LinearInterpolations.interpolate(modeDampings[c,i,mode][jO:jO+1],URange[jO:jO+1],0))
                    push!(flutterOnsetFreqOfMode[c,i,mode],LinearInterpolations.interpolate(modeDampings[c,i,mode][jO:jO+1],modeFrequencies[c,i,mode][jO:jO+1],0))
                    push!(flutterOnsetTipOOPOfMode[c,i,mode], tipOOP[c,i,jO])
                end
            end
        end
        # All flutter onset speeds and frequencies and tip OOP
        flutterOnsetSpeedsAll[c,i] = vcat(filter(!isempty,flutterOnsetSpeedOfMode[c,i,:])...)
        flutterOnsetFreqsAll[c,i] = vcat(filter(!isempty,flutterOnsetFreqOfMode[c,i,:])...)
        flutterOnsetTipOOPAll = vcat(filter(!isempty,flutterOnsetTipOOPOfMode[c,i,:])...)
        # Lowest flutter onset speed, corresponding frequency and tip OOP
        if !isempty(flutterOnsetSpeedsAll[c,i])
            iLowest = sortperm(flutterOnsetSpeedsAll[c,i])[1]
            flutterOnsetSpeed[c,i] = flutterOnsetSpeedsAll[c,i][iLowest]
            flutterOnsetFreq[c,i] = flutterOnsetFreqsAll[c,i][iLowest]
            flutterOnsetTipOOP[c,i] = flutterOnsetTipOOPAll[iLowest]
        end
    end
end

# Set paths
relPath = "/dev/sweptPazy/S10/outputs/S10PazyFlutterOnsetPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
configColors = cgrad(:rainbow, length(tipMassConfigs), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 1
gr()

# Flutter onset speed vs. tip displacement
plt = plot(xlabel="Tip displacement [m]", ylabel="Flutter speed [m/s]", xlims=[0,0.3+0.01], ylims=[0,80], xticks=vcat(0:0.05:0.3), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
# plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
plot!([NaN], [NaN], lw=lw, c=configColors[1], label="LE weight")
plot!([NaN], [NaN], lw=lw, c=configColors[2], label="TE weight")
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:], flutterOnsetSpeed[c,:], c=configColors[c], lw=lw, label=false)
end
display(plt)
savefig(plt,string(absPath,"/S10PazyFlutterOnsetPitchRange_Uf_vs_tipDisp.pdf"))

println("Finished S10PazyFlutterPitchRange.jl")