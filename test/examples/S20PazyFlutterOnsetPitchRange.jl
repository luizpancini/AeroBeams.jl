using AeroBeams, DelimitedFiles, LinearInterpolations

# Tip mass configurations (assumed tipMassPosOffsetLE/tipMassPosOffsetTE ahead/behind the LE/TE)
tipMassConfigs = ["LE", "TE"]
tipMass = 15.5e-3
tipMassPosOffsetLE = 0.05
tipMassPosOffsetTE = 0.05

# Root pitch angle range
θRange = π/180*vcat(0.5,1:1:12)

# Airspeed range
URange = collect(1:1:85)

# Sweep angle [rad]
Λ = 20*π/180

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Gravity
g = 9.80665

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

# Sweep tip weight configurations
for (c,config) in enumerate(tipMassConfigs)
    tipMassPos = config == "LE" ? (chord*normSparPos + tipMassPosOffsetLE) : (- (chord*(1-normSparPos) + tipMassPosOffsetTE))
    # Sweep root pitch angle
    for (i,θ) in enumerate(θRange)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for tip mass at $config, θ = $(round(θ*180/pi,digits=1)) deg, U = $U m/s")
            # Model
            model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,g=g,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0])
            # Create and solve problem
            problem = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[2π,Inf])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[c,i,j] = problem.frequenciesOscillatory
            untrackedDamps[c,i,j] = problem.dampingsOscillatory
            untrackedEigenvectors[c,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP deflection
            tipOOP[c,i,j] = problem.nodalStatesOverσ[end][nElem].u_n2_b[3]
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
                    push!(flutterOnsetSpeedOfMode[c,i,mode],interpolate(modeDampings[c,i,mode][jO:jO+1],URange[jO:jO+1],0))
                    push!(flutterOnsetFreqOfMode[c,i,mode],interpolate(modeDampings[c,i,mode][jO:jO+1],modeFrequencies[c,i,mode][jO:jO+1],0))
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

# Load reference data (experiments from Technion - https://doi.org/10.5281/zenodo.16354530)
Uf_vs_tipdisp_LE_c = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_LE_c.txt")
Uf_vs_tipdisp_LE_g = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_LE_g.txt")
Uf_vs_tipdisp_TE_c = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_TE_c.txt")
Uf_vs_tipdisp_TE_g = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_TE_g.txt")

println("Finished S20PazyFlutterPitchRange.jl")