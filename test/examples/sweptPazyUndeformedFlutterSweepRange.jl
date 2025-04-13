using AeroBeams

# Sweep angle range [rad]
ΛRange = [0,10,20,30]*π/180

# Airspeed range
URange = collect(1:1:120)

# Root pitch angle
θ = 0*π/180

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

# Number of modes
nModes = 3

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter,alwaysUpdateJacobian=false)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(ΛRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
tipOOP = Array{Float64}(undef,length(ΛRange),length(URange))

# Sweep angle of sweep
for (i,Λ) in enumerate(ΛRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for Λ = $(round(Λ*180/π)) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,g=g,sweepStructuralCorrections=sweepStructuralCorrections)
        # Create and solve problem
        lowerLimit = U <= 80 ? 2π : min(2π*U/20,2π*4)
        problem = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[lowerLimit,Inf],systemSolver=NR)
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

# Load reference data (from AePW4 meetings)
undef_damp_Lambda0_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_damp_Lambda0_NastranTechnion.txt")
undef_damp_Lambda20_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_damp_Lambda20_NastranTechnion.txt")
undef_damp_Lambda30_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_damp_Lambda30_NastranTechnion.txt")
undef_freq_Lambda0_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_freq_Lambda0_NastranTechnion.txt")
undef_freq_Lambda20_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_freq_Lambda20_NastranTechnion.txt")
undef_freq_Lambda30_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_freq_Lambda30_NastranTechnion.txt")

println("Finished sweptPazyUndeformedFlutterSweepRange.jl")