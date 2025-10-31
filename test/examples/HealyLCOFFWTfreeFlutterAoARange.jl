using AeroBeams, DelimitedFiles

# Wingtip configuration
wingtipTab = false

# Pitch angle range [rad]
θRange = π/180*[2.5; 7.5]

# Airspeed range [m/s]
URange = collect(5:0.5:30)

# Stiffness of the spring around the hinge
kSpring = 1e-4
kIPBendingHinge = 1e4

# Discretization
nElementsInner = 8
nElementsOuter = 3
nElementsFFWT = 4

# Tip loss flag
hasTipCorrection = true
tipLossDecayFactor = nothing

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Number of modes
nModes = 10

# Initialize outputs
untrackedFreqs = [fill(NaN64, nModes) for _ in θRange, _ in URange]
untrackedDamps = [fill(NaN64, nModes) for _ in θRange, _ in URange]
untrackedEigenvectors = [fill(NaN64 + im*NaN64, nModes, nModes) for _ in θRange, _ in URange]
freqs = [fill(NaN64, nModes) for _ in θRange, _ in URange]
damps = [fill(NaN64, nModes) for _ in θRange, _ in URange]
modeFrequencies = [fill(NaN64, nModes) for _ in θRange, _ in URange]
modeDampings = [fill(NaN64, nModes) for _ in θRange, _ in URange]
modeDampingRatios = [fill(NaN64, nModes) for _ in θRange, _ in URange]
ϕHinge = fill(NaN64, length(θRange), length(URange))
problem = Array{EigenProblem}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (j,θ) in enumerate(θRange)
    # Sweep airspeed
    for (k,U) in enumerate(URange)
        # Display progress 
        println("Solving for θ=$(round(θ*180/π,digits=1)) deg, U=$U m/s") 
        # Update model
        model = create_HealyLCOFFWT(wingtipTab=wingtipTab,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)
        # Set initial guess solution as the one from previous airspeed
        x0 = (k==1 || !problem[j,k-1].systemSolver.convergedFinalSolution) ? zeros(0) : problem[j,k-1].x
        # Create and solve problem
        problem[j,k] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,frequencyFilterLimits=[1e-2*U,Inf],x0=x0)
        solve!(problem[j,k])
        # Get outputs, if converged
        if problem[j,k].systemSolver.convergedFinalSolution
            # Frequencies and dampings
            untrackedFreqs[j,k] = problem[j,k].frequenciesOscillatory
            untrackedDamps[j,k] = problem[j,k].dampingsOscillatory
            untrackedEigenvectors[j,k] = problem[j,k].eigenvectorsOscillatoryCplx
            # Hinge fold angle
            ϕHinge[j,k] = -problem[j,k].model.hingeAxisConstraints[1].ϕ*180/π
        end
    end
    # Apply mode tracking
    freqs[j,:],damps[j,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[j,:],untrackedDamps[j,:],untrackedEigenvectors[j,:])
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[j,mode] = [freqs[j,k][mode] for k in eachindex(URange)]
        modeDampings[j,mode] = [damps[j,k][mode] for k in eachindex(URange)]
        modeDampingRatios[j,mode] = modeDampings[j,mode]./modeFrequencies[j,mode]
    end
end

println("Finished HealyBaselineFFWTfreeFlutterAoARangeURange.jl")