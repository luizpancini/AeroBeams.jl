using AeroBeams

# Airspeed range [m/s]
URange = collect(10:1:40)

# Flare angle [rad]
Λ = 30*π/180

# Pitch angle [rad]
θ = 3*π/180

# Sideslip angle [rad]
β = 0*π/180

# Spring stiffness [Nm/rad]
kSpring = 1e-4

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution, especially at lower airspeeds)
withTipCorrection = true
tipLossDecayFactor = 10

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-6
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

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
    model = create_HealyFFWT(flareAngle=Λ,kSpring=kSpring,airspeed=U,pitchAngle=θ,withTipCorrection=withTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,flightDirection=[sin(β);cos(β);0])
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
freqs,damps,_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished HealyFFWTflutterURange.jl")