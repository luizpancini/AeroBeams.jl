using AeroBeams

# Sideslip angle range [rad]
βRange = π/180*vcat(-20:1:30)

# Flare angle [rad]
Λ = 30*π/180

# Pitch angle [rad]
θ = 3*π/180

# Airspeed [m/s]
U = 33

# Spring stiffness [Nm/rad]
kSpring = 1e-4

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution, especially at lower airspeeds)
hasTipCorrection = true
tipLossDecayFactor = 10

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Number of modes
nModes = 3

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(βRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(βRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(βRange))
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
problem = Array{EigenProblem}(undef,length(βRange))

# Sweep sideslip angle
for (i,β) in enumerate(βRange)
    # Display progress
    println("Solving for β=$(round(Int,β*180/π)) deg")
    # Update model
    model = create_HealyFFWT(flareAngle=Λ,kSpring=kSpring,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,flightDirection=[sin(β);cos(β);0])
    # Set initial guess solution as the one from previous sideslip angle
    x0 = (i>1 && problem[i-1].systemSolver.convergedFinalSolution) ? problem[i-1].x : zeros(0)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,x0=x0,frequencyFilterLimits=[1e0,Inf])
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(βRange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(βRange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(βRange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished HealyFFWTflutterSideslipRange.jl")