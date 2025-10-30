using AeroBeams

# Airspeed range [m/s]
URange = collect(1:1:40)

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

# Spring stiffness [Nm/rad]
kSpring = 1e-4
kIPBendingHinge = 1e1

# Solution method for hinge constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = true

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
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
    model = create_PazyFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,hingeNode=hingeNode,flareAngle=Λ,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,foldAngle=foldAngle,flightDirection=[sin(β);cos(β);0],tipMass=tipMass,ηtipMass=ηtipMass)
    # Initial guess solution
    x0 = i==1 ? zeros(0) : problem[i-1].x
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,frequencyFilterLimits=[1,Inf],x0=x0)
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = problem[i].dampingsOscillatory
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
    # Coast (or fold) angle
    ϕHinge[i] = problem[i].model.hingeAxisConstraints[1].ϕ*180/π
end

# Apply mode tracking
freqs,damps,_ = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished PazyFFWTflutterURange.jl")