using AeroBeams, LinearInterpolations

# Option for mode tracking
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Derivation method
derivationMethod = AD()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Set system solver options
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of modes
nModes = 5

# Set pitch angle and airspeed ranges, and initialize outputs
θRange = collect(-0.25:0.25:7)
URange = collect(0:0.5:90)
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
flutterOnsetSpeedsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOnsetFreqsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOnsetDispOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOffsetSpeedsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOffsetFreqsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOffsetDispOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Model
        PazyWingFlutterPitchRange,nElem,L,chord,normSparPos = create_Pazy(aeroSolver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,upright=upright,θ=θ*π/180,airspeed=U)
        # Create and solve problem
        problem = create_EigenProblem(model=PazyWingFlutterPitchRange,nModes=nModes,systemSolver=NR)
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem.dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        # Get OOP displacement at midchord
        tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist = asind(Δ[3])
        tip_OOP[i,j] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist))
    end
    # Apply mode tracking, if applicable
    if modeTracking
        freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and damping ratios by mode
    modeFrequencies = Array{Vector{Float64}}(undef,nModes)
    modeDampings = Array{Vector{Float64}}(undef,nModes)
    modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
    for mode in 1:nModes
        modeFrequencies[mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
    end
    # Loop over modes: compute flutter onset and offset speeds, respective frequencies and OOP displacements
    for mode in 1:nModes
        # Find flutter onset indices
        onsetIndices = findall((modeDampings[mode][2:end] .> 0) .& (modeDampings[mode][1:end-1] .< 0)) .+ 1 
        nIndOn = length(onsetIndices)
        # Loop flutter onset indices
        flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn)
        for (n,k) in enumerate(onsetIndices)
            flutterOnsetSpeeds[n] = interpolate(modeDampings[mode][k-1:k],URange[k-1:k],0)
            flutterOnsetFreqs[n] = interpolate(modeDampings[mode][k-1:k],modeFrequencies[mode][k-1:k],0)
            flutterOnsetDisp[n] = interpolate(modeDampings[mode][k-1:k],tip_OOP[i,k-1:k]/L*100,0)
        end
        if nIndOn == 0
            flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = [NaN],[NaN],[NaN]
        end
        # Set flutter onset variables for current mode
        flutterOnsetSpeedsOfMode[i,mode] = flutterOnsetSpeeds
        flutterOnsetFreqsOfMode[i,mode] = flutterOnsetFreqs
        flutterOnsetDispOfMode[i,mode] = flutterOnsetDisp
        # Find flutter offset indices
        offsetIndices = findall((modeDampings[mode][2:end] .< 0) .& (modeDampings[mode][1:end-1] .> 0)) .+ 1
        nIndOff = length(offsetIndices)
        # Find flutter offset variables
        flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff)
        # Loop flutter offset indices
        for (n,k) in enumerate(offsetIndices)
            flutterOffsetSpeeds[n] = interpolate(-modeDampings[mode][k-1:k],URange[k-1:k],0)
            flutterOffsetFreqs[n] = interpolate(-modeDampings[mode][k-1:k],modeFrequencies[mode][k-1:k],0)
            flutterOffsetDisp[n] = interpolate(-modeDampings[mode][k-1:k],tip_OOP[i,k-1:k]/L*100,0)
        end
        if nIndOff == 0
            flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = [NaN],[NaN],[NaN]
        end
        # Set flutter offset variables for current mode
        flutterOffsetSpeedsOfMode[i,mode] = flutterOffsetSpeeds
        flutterOffsetFreqsOfMode[i,mode] = flutterOffsetFreqs
        flutterOffsetDispOfMode[i,mode] = flutterOffsetDisp
    end
end

# Load reference data from DRACHINSKI et al. "Flutter Tests of the Pazy Wing"(2021)
rootPitchVelUp = [3; 5; 7]
rootPitchVelDown = [3; 5]
flutterOnsetVelUp = [49; 43; 38]
flutterOffsetVelUp = [58; 51; 46]
flutterOnsetVelDown = [55; 48]
flutterOffsetVelDown = [40; 36]
flutterOnsetDispUp = [23; 24.5; 26]
flutterOnsetDispDown = [28.5; 31.5]

println("Finished PazyWingFlutterPitchRange.jl")