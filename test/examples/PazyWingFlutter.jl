using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = BLi()

# Derivation method
derivationMethod = AD()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Root pitch angle [rad]
θ = 5*π/180

# Option for mode tracking
modeTracking = true

# Initialize model
PazyWingFlutter,nElem,L,chord,normSparPos = create_Pazy(aeroSolver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,upright=upright,θ=θ)

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,alwaysUpdateJacobian=false)

# Number of modes
nModes = 5

# Set airspeed range, and initialize outputs
URange = collect(0:0.5:100)
problem = Array{EigenProblem}(undef,length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
tip_OOP = Array{Float64}(undef,length(URange))
u1_of_x1 = Array{Vector{Float64}}(undef,length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))

# Undeformed nodal and midpoint positions
x1_0 = vcat([vcat(PazyWingFlutter.beams[1].elements[e].r_n1[3],PazyWingFlutter.beams[1].elements[e].r_n2[3]) for e in 1:nElem]...)
x3_0 = -vcat([vcat(PazyWingFlutter.beams[1].elements[e].r_n1[1],PazyWingFlutter.beams[1].elements[e].r_n2[1]) for e in 1:nElem]...)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=PazyWingFlutter,v_A=[0;U;0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=PazyWingFlutter,nModes=nModes,systemSolver=NR,frequencyFilterLimits=[1.5,Inf64])
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
    # Get OOP displacement at midchord
    tip_p = problem[i].nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist = asind(Δ[3])
    tip_OOP[i] = -(problem[i].nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist))
    # Displacements over span
    u1_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1_b[1],problem[i].nodalStatesOverσ[end][e].u_n2_b[1]) for e in 1:nElem]...)
    u3_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1_b[3],problem[i].nodalStatesOverσ[end][e].u_n2_b[3]) for e in 1:nElem]...)
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1_of_x1[i]
    x3_def[i] = x3_0 .+ u3_of_x1[i]
end

# Apply mode tracking, if applicable
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Update frequencies and dampings order on problem
for i in eachindex(URange)
    problem[i].frequenciesOscillatory = freqs[i]
    problem[i].dampingsOscillatory = damps[i]
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

println("Finished PazyWingFlutter.jl")