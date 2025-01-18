using AeroBeams

# Hinge configuration
hingeConfiguration = "locked"

# Airspeed range [m/s]
URange = collect(0:1:40)

# Pitch angle [rad]
θ = 0*π/180

# Discretization
nElementsInner = 16
nElementsFFWT = 4

# Tip loss options (assumed, since Healy's analysis uses DLM for aerodynamic)
withTipCorrection = true
tipLossDecayFactor = 12

# System solver
σ0 = 1
maxIter = 50
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Number of modes
nModes = 4

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
problem = Array{EigenProblem}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U=$U m/s")
    # Update model
    model =  create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,airspeed=U,pitchAngle=θ,withTipCorrection=withTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,frequencyFilterLimits=[1e-2,Inf])
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished HealyBaselineFFWTlockedFlutterURange.jl")