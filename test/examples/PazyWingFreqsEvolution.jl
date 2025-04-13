using AeroBeams, DelimitedFiles

# Aerodynamic solver
aeroSolver = Indicial()

# Flag for skin (wing covering)
withSkin = false

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "Exponential"

# Flag for upright position
upright = true

# Airspeed range
URange = collect(0:0.5:60)

# Pitch angle [rad]
θ = 7*π/180

# Gravity
g = 0

# Flag for small angles approximation
smallAngles = true

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Set system solver options
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of vibration modes
nModes = 5

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
tipOOP = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model
    PazyWingFreqsEvolution,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,withSkin=withSkin,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,upright=upright,θ=θ,airspeed=U,g=g,smallAngles=smallAngles)
    # Create and solve problem
    problem = create_EigenProblem(model=PazyWingFreqsEvolution,nModes=nModes,systemSolver=NR)
    solve!(problem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem.dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem.eigenvectorsOscillatoryCplx
    # Get OOP displacement at midchord
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist = asind(Δ[3])
    tipOOP[i] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist))
end

# Apply mode tracking
freqs,_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end

# Load reference data
freqsVsU_alpha7 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/freqsVsU_alpha7.txt")
freqsVsDisp_alpha7 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/freqsVsDisp_alpha7.txt")

println("Finished PazyWingFreqsEvolution.jl")