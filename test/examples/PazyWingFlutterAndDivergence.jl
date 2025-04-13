using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Inflow()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Gravity
g = 0

# Pitch angle
θ = 0*π/180

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass
tipMass = 0.00
ηtipMass = [0; -chord*(1-normSparPos)-0.00; 0]

# Set system solver options
σ0 = 1
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of modes
nModes = 3

# Set airspeed range, and initialize outputs
URange = collect(5:0.5:120)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
nonOscillatoryDampings = Array{Vector{Float64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
tip_OOP = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model
    PazyWingFlutterAndDivergence,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,θ=θ,airspeed=U,g=g,tipMass=tipMass,ηtipMass=ηtipMass)
    # Create and solve problem
    problem = create_EigenProblem(model=PazyWingFlutterAndDivergence,nModes=nModes,systemSolver=NR)
    solve!(problem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem.dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem.eigenvectorsOscillatoryCplx
    nonOscillatoryDampings[i] = problem.dampingsNonOscillatory
    # Get OOP displacement at midchord
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist = asind(Δ[3])
    tip_OOP[i] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist))
end

# Apply mode tracking, if applicable
freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

# Flutter onset and offset speeds, respective frequencies and OOP displacements
flutterOnsetSpeedsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOnsetFreqsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOnsetDispOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOffsetSpeedsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOffsetFreqsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOffsetDispOfMode = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    # Find flutter onset indices
    onsetIndices = findall((modeDampings[mode][2:end] .> 0) .& (modeDampings[mode][1:end-1] .< 0)) .+ 1 
    nIndOn = length(onsetIndices)
    # Loop flutter onset indices
    flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn)
    for (n,i) in enumerate(onsetIndices)
        flutterOnsetSpeeds[n] = interpolate(modeDampings[mode][i-1:i],URange[i-1:i],0)
        flutterOnsetFreqs[n] = interpolate(modeDampings[mode][i-1:i],modeFrequencies[mode][i-1:i],0)
        flutterOnsetDisp[n] = interpolate(modeDampings[mode][i-1:i],tip_OOP[i-1:i]/L*100,0)
        println("Flutter onset speed = $(flutterOnsetSpeeds[n]) m/s, flutter frequency = $(flutterOnsetFreqs[n]/(2*π)) Hz, flutter OOP disp = $(flutterOnsetDisp[n]) % semispan (mode $mode)")
    end
    # Set flutter onset variables for current mode
    flutterOnsetSpeedsOfMode[mode] = flutterOnsetSpeeds
    flutterOnsetFreqsOfMode[mode] = flutterOnsetFreqs
    flutterOnsetDispOfMode[mode] = flutterOnsetDisp
    # Find flutter offset indices
    offsetIndices = findall((modeDampings[mode][2:end] .< 0) .& (modeDampings[mode][1:end-1] .> 0)) .+ 1
    nIndOff = length(offsetIndices)
    # Find flutter offset variables
    flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff)
    # Loop flutter offset indices
    for (n,i) in enumerate(offsetIndices)
        flutterOffsetSpeeds[n] = interpolate(-modeDampings[mode][i-1:i],URange[i-1:i],0)
        flutterOffsetFreqs[n] = interpolate(-modeDampings[mode][i-1:i],modeFrequencies[mode][i-1:i],0)
        flutterOffsetDisp[n] = interpolate(-modeDampings[mode][i-1:i],tip_OOP[i-1:i]/L*100,0)
        println("Flutter offset speed = $(flutterOffsetSpeeds[n]) m/s, flutter frequency = $(flutterOffsetFreqs[n]/(2*π)) Hz, flutter OOP disp = $(flutterOffsetDisp[n]) % semispan (mode $mode)")
    end
    # Set flutter offset variables for current mode
    flutterOffsetSpeedsOfMode[mode] = flutterOffsetSpeeds
    flutterOffsetFreqsOfMode[mode] = flutterOffsetFreqs
    flutterOffsetDispOfMode[mode] = flutterOffsetDisp
end

# Initialize divergence speed, its index, and flag for divergence being found
global divergenceFound = nothing
global iD = nothing
global divergenceFound = false

# Separate non-oscillatory dampings by mode
nNOModes = 1
modeNonOscillatoryDampings = Array{Vector{Float64}}(undef,nNOModes)
modeNonOscillatoryDampingsEst = Array{Vector{Float64}}(undef,nNOModes)
for mode in 1:nNOModes
    # Mode dampings
    modeNonOscillatoryDampings[mode] = [nonOscillatoryDampings[i][mode] for i in eachindex(URange)]
    # Estimated mode dampings from backward finite difference extrapolation
    modeNonOscillatoryDampingsEst[mode] = backward_extrapolation(modeNonOscillatoryDampings[mode])
    # Divergence is found when the sign of the estimated value is different from the actual
    if !divergenceFound
        global iD = findfirst(i -> modeNonOscillatoryDampings[mode][i]*modeNonOscillatoryDampingsEst[mode][i] < 0, 1:length(URange))
        if !isnothing(iD)
            global divergenceSpeed = LinearInterpolations.interpolate(modeNonOscillatoryDampingsEst[mode][iD-1:iD],URange[iD-1:iD],0)
            global divergenceFound = true
            println("Divergence speed = $divergenceSpeed m/s")
        end
    end
end

println("Finished PazyWingFlutterAndDivergence.jl")