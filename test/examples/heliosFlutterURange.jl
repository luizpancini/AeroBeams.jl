using AeroBeams

# Stiffness factor
λ = 1

# TF to include beam pods
beamPods = false

# Aerodynamic solver
aeroSolver = Indicial()

# Payload [lb]
P = 400

# Option for mode tracking
modeTracking = true

# Number of modes
nModes = 10

# Set NR system solver for trim problem
relaxFactor = 0.5
σ0 = 1
σstep = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter,displayStatus=false)

# Set airspeed range, and initialize outputs
URange = collect(10.5:0.5:30)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model for trim problem
    heliosTrim,midSpanElem,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
    # Add springs at wing root
    add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
    # Set initial guess solution as previous known solution
    x0Trim = (i==1) ? zeros(0) : trimProblem.x
    # Create and solve trim problem
    global trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
    trimδ = trimProblem.x[end]
    # Model for eigen problem
    heliosEigen,_,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust)
    # Add springs at wing root
    add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
    # Set initial solution as trim solution
    x0Eigen = trimProblem.x[1:end-2]
    # Create and solve eigen problem
    global eigenProblem = create_EigenProblem(model=heliosEigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],jacobian=trimProblem.jacobian[1:end,1:end-2],inertia=trimProblem.inertia)
    solve_eigen!(eigenProblem)
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

println("Finished heliosFlutterURange.jl")