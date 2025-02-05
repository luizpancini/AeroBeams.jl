using AeroBeams

# Airspeed range [m/s]
URange = collect(1:1:60)

# Hinge node
hingeNode = 13

# Fold angle [rad]
foldAngle = nothing

# Flare angle [rad]
Λ = 10*π/180

# Pitch angle [rad]
θ = 3*π/180

# Sideslip angle [rad]
β = 0*π/180

# Tip mass [kg] and its position [m]
tipMass = 10e-3
ηtipMass = 1e-3*[0;50;0]

# Spring stiffness [Nm/rad] (as kSpring → ∞, the hingeless behavior is obtained, as expected!)
kSpring = 1e-4

# Solution method for hinge constraint
solutionMethod = "appliedMoment"

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-6
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Number of modes
nModes = 5

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
ϕHinge = Array{Float64}(undef,length(URange))
problem = Array{EigenProblem}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U=$U m/s")
    # Update model
    model = create_PazyFFWT(hingeNode=hingeNode,flareAngle=Λ,kSpring=kSpring,airspeed=U,pitchAngle=θ,foldAngle=foldAngle,flightDirection=[sin(β);cos(β);0],tipMass=tipMass,tipMassPosition=ηtipMass)
    # Set initial guess solution as the one from previous sideslip angle
    x0 = (i>1 && problem[i-1].systemSolver.convergedFinalSolution) ? problem[i-1].x : zeros(0)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,x0=x0,frequencyFilterLimits=[1e-0,Inf])
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
    # Coast (or fold) angle
    ϕHinge[i] = problem[i].model.hingeAxisConstraints[1].ϕ*180/π
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished PazyFFWTflutterURange.jl")