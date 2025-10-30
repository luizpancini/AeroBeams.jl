using AeroBeams

# Sideslip angle range [rad]
βRange = π/180*vcat(-20:1:20)

# Pitch angle range [rad]
θRange = π/180*[-3,3,9,15]

# Flare angle [rad]
Λ = 30*π/180

# Airspeed [m/s]
U = 25

# Solution method for hinge constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = true

# Spring stiffness [Nm/rad]
kSpring = 1e-4
kIPBendingHinge = 1e-1

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 12

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Number of modes
nModes = 2

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(βRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(βRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(βRange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(βRange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(βRange))
modeFrequencies = Array{Vector{Float64}}(undef,length(θRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(θRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(θRange),nModes)
problem = Array{EigenProblem}(undef,length(θRange),length(βRange))

# Sweep root pitch angle
for (i,θ) in enumerate(θRange)
    # Sweep sideslip angle
    for (j,β) in enumerate(βRange)
        # Display progress
        println("Solving for θ=$(round(θ*180/π,digits=1)) deg, β=$(round(β*180/π,digits=1)) deg")
        # Update model
        model = create_HealySideslipFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,flareAngle=Λ,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,flightDirection=[-sin(β);cos(β);0])
        # Set initial guess solution as the one from previous sideslip angle
        x0 = (j>1 && problem[i,j-1].systemSolver.convergedFinalSolution) ? problem[i,j-1].x : zeros(0)
        # Create and solve problem
        problem[i,j] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,x0=x0,frequencyFilterLimits=[1e0,Inf])
        solve!(problem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = problem[i,j].dampingsOscillatory
        untrackedEigenvectors[i,j] = problem[i,j].eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(βRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,k][mode] for k in eachindex(βRange)]
        modeDampings[i,mode] = [damps[i,k][mode] for k in eachindex(βRange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

println("Finished HealySideslipFFWTflutterSideslipRange.jl")