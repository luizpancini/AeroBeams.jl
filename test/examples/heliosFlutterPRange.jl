using AeroBeams

# Wing airfoil
wingAirfoil = NACA23012A

# Option for reduced chord
reducedChord = false

# Option to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Stiffness factor
λ = 1

# Aerodynamic solver
aeroSolver = BLi()

# Airspeed
U = 40*0.3048

# Number of modes
nModes = 10

# System solver for trim problem
relaxFactor = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=100,pseudoInverseMethod=:MoorePenrose,displayStatus=false)

# Set payload range, and initialize outputs
PRange = collect(0:20:500)
untrackedFreqs = Array{Vector{Float64}}(undef,length(PRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(PRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(PRange))
freqs = Array{Vector{Float64}}(undef,length(PRange))
damps = Array{Vector{Float64}}(undef,length(PRange))

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)

# Sweep payload
for (i,P) in enumerate(PRange)
    # Display progress
    println("Solving for payload = $P lb")
    # Model for trim problem
    heliosTrim,midSpanElem,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,payloadOnWing=payloadOnWing)
    # Add springs at wing root
    add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
    # Update model
    update_model!(heliosTrim)
    # Set initial guess solution as previous known solution
    x0Trim = (i==1) ? zeros(0) : trimProblem.x
    # Create and solve trim problem
    global trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimAoA = trimProblem.aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
    trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
    trimδ = trimProblem.x[end]
    println("AoA = $(trimAoA), T = $(trimThrust), δ = $(trimδ*180/π)")
    # Model for eigen problem
    heliosEigen,_,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,reducedChord=reducedChord,payloadOnWing=payloadOnWing)
    # Add springs at wing root
    add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
    # Update model
    update_model!(heliosEigen)
    # Create and solve eigen problem
    global eigenProblem = create_EigenProblem(model=heliosEigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],refTrimProblem=trimProblem)
    solve_eigen!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = eigenProblem.dampingsOscillatory
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
end

# Apply mode tracking
freqs,damps,_,matchedModes = mode_tracking_hungarian(PRange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(PRange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(PRange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

println("Finished heliosFlutterPRange.jl")