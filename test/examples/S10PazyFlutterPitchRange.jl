using AeroBeams, DelimitedFiles

# Tip mass configurations (assumed tipMassPosOffsetLE/tipMassPosOffsetTE ahead/behind the LE/TE)
tipMassConfigs = ["LE", "TE"]
tipMass = 11e-3
tipMassPosOffsetLE = 0.01
tipMassPosOffsetTE = 0.05

# Root pitch angle range
θRange = π/180*[0,1,3,5,7,10]

# Airspeed range
URange = collect(1:1:100)

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
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Gravity
g = 9.80665

# Number of modes
nModes = 5

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

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
    # Sweep pitch angle
    for (i,θ) in enumerate(θRange)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for for tip mass at $config, θ = $(round(θ*180/π)) deg, U = $U m/s")
            # Model
            sweptPazyFlutterPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,g=g,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0])
            # Create and solve problem
            problem = create_EigenProblem(model=sweptPazyFlutterPitchRange,nModes=nModes,frequencyFilterLimits=[2π,Inf])
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

println("Finished sweptPazyFlutterPitchRange.jl")