using AeroBeams

# Stiffness factor
λ = 1

# Option for mode tracking
modeTracking = true

# TF to include beam pod
beamPods = false

# Aerodynamic solver
aeroSolver = Indicial()

# Number of modes
nModes = 8

# Set system solver options
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set airspeed range, and initialize outputs
URange = collect(0.1:0.1:15)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model 
    _,_,heliosWing,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,airspeed=U)
    # Create and solve trim problem
    global eigenProblem = create_EigenProblem(model=heliosWing,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],systemSolver=NR)
    solve!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(eigenProblem.dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
end

# Apply mode tracking, if applicable
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished heliosWingFlutter.jl")