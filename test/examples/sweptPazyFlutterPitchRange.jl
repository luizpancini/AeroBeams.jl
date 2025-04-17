using AeroBeams

# Sweep angle [rad]
Λ = 20*π/180

# Root pitch angle range
θRange = π/180*[0,1,3,5,7]

# Airspeed range
URange = collect(0.5:0.5:100)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = false

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass and its position over the chord
tipMass = 0.0
tipMassPos = 0 # EA
# tipMassPos = chord*normSparPos # LE
# tipMassPos = -chord*(1-normSparPos) # TE

# Number of modes
nModes = 3

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter,alwaysUpdateJacobian=false)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(θRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(θRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(θRange),nModes)
tipOOP = Array{Float64}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(θ*180/π)) deg, U = $U m/s")
        # Model
        sweptPazyFlutterPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,g=g,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0])
        # Create and solve problem
        problem = create_EigenProblem(model=sweptPazyFlutterPitchRange,nModes=nModes,frequencyFilterLimits=[2π,Inf],systemSolver=NR)
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem.dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        # Tip OOP deflection
        tipOOP[i,j] = problem.nodalStatesOverσ[end][nElem].u_n2_b[3]
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

println("Finished sweptPazyFlutterPitchRange.jl")