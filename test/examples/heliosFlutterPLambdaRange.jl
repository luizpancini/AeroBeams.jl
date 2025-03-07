using AeroBeams

# Wing airfoil
wingAirfoil = HeliosWingAirfoil

# Option for reduced chord
reducedChord = false

# Option to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Airspeed
U = 40*0.3048

# Option for mode tracking
modeTracking = true

# Number of modes
nModes = 10

# System solver for trim problem
relaxFactor = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=100,displayStatus=false)

# Set stiffness factor and payload ranges, and initialize outputs
λRange = [1,5,100]
PRange = collect(0:10:500)
untrackedFreqs = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(λRange),length(PRange))
freqs = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
damps = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
modeFrequencies =  Array{Vector{Float64}}(undef,length(λRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(λRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(λRange),nModes)
eigenProblem = Array{EigenProblem}(undef,length(λRange),length(PRange))

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Solving for λ = $λ, payload = $P lb")
        # Model for trim problem
        heliosTrim,midSpanElem,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,payloadOnWing=payloadOnWing)
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        heliosTrim.skipValidationMotionBasisA = true
        update_model!(heliosTrim)
        # Set initial guess solution as previous known solution
        x0Trim = (j==1) ? zeros(0) : trimProblem.x
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
        heliosEigen.skipValidationMotionBasisA = true
        update_model!(heliosEigen)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=heliosEigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],jacobian=trimProblem.jacobian[1:end,1:end-2],inertia=trimProblem.inertia,refTrimProblem=trimProblem)
        solve_eigen!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-12)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking, if applicable
    if modeTracking
        freqs[i,:],damps[i,:],_,matchedModes = mode_tracking(PRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(PRange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(PRange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

println("Finished heliosFlutterPLambdaRange.jl")