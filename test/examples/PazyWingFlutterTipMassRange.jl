using AeroBeams

# Option for mode tracking
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Derivation method
derivationMethod = AD()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Pitch angle
θ = 0*π/180

# Set system solver options
σ0 = 1.0
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of modes
nModes = 3

# Get wing chord and spar position
_,_,chord,normSparPos = geometrical_properties_Pazy()

# Set tip mass and its position ranges, airspeed range, and initialize outputs
configurations = [1; 2; 3]
tipMassRange = [0; 0.01; 0.01]
tipMassPosRange = [0; chord*(1-normSparPos); -chord*normSparPos]
URange = collect(0:0.5:100)
untrackedFreqs = Array{Vector{Float64}}(undef,3,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,3,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,3,length(URange))
freqs = Array{Vector{Float64}}(undef,3,length(URange))
damps = Array{Vector{Float64}}(undef,3,length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,3,nModes)
modeDampings = Array{Vector{Float64}}(undef,3,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,3,nModes)

# Sweep configurations
for c in configurations   
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for configuration $c, U = $U m/s")
        # Model
        PazyWingFlutterTipMassRange,_ = create_Pazy(aeroSolver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,upright=upright,θ=θ,airspeed=U,tipMass=tipMassRange[c],ξtipMass=[0;tipMassPosRange[c];0])
        # Create and solve problem
        problem = create_EigenProblem(model=PazyWingFlutterTipMassRange,nModes=nModes,systemSolver=NR)
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[c,j] = problem.frequenciesOscillatory
        untrackedDamps[c,j] = round_off!(problem.dampingsOscillatory,1e-8)
        untrackedEigenvectors[c,j] = problem.eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking, if applicable
    if modeTracking
        freqs[c,:],damps[c,:],_ = mode_tracking(URange,untrackedFreqs[c,:],untrackedDamps[c,:],untrackedEigenvectors[c,:])
    else
        freqs[c,:],damps[c,:] = untrackedFreqs[c,:],untrackedDamps[c,:]
    end
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[c,mode] = [freqs[c,i][mode] for i in eachindex(URange)]
        modeDampings[c,mode] = [damps[c,i][mode] for i in eachindex(URange)]
        modeDampingRatios[c,mode] = modeDampings[c,mode]./modeFrequencies[c,mode]
    end
end

println("Finished PazyWingFlutterTipMassRange.jl")