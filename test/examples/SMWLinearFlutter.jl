using AeroBeams, LinearInterpolations

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Altitude
h = 20e3

# Gravity
g = 9.80665

# Pitch angle
θ = 0

# Discretization
nElem = 16

# Model
SMWLinearFlutter,L = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ,nElem=nElem,altitude=h,g=g)

# Set airspeed range and initialize outputs
URange = collect(0:0.1:40)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
problem = Array{EigenProblem}(undef,length(URange))

# Set number of vibration modes
nModes = 5

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update velocity of basis A 
    set_motion_basis_A!(model=SMWLinearFlutter,v_A=[0;U;0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=SMWLinearFlutter,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],getLinearSolution=true)
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
end

# Frequencies and dampings after mode tracking
freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

# Find flutter speed and flutter frequency 
global flutterSpeed,flutterFreq = NaN,NaN
dampsOfMode = Array{Vector{Float64}}(undef,nModes)
freqsOfMode = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    dampsOfMode[mode] = [damps[j][mode] for j in eachindex(URange)]
    freqsOfMode[mode] = [freqs[j][mode] for j in eachindex(URange)]
    indexInstability = findfirst(x->x>0,dampsOfMode[mode])
    if isnothing(indexInstability)
        continue
    end
    global flutterSpeed = interpolate(dampsOfMode[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
    global flutterFreq = interpolate(dampsOfMode[mode][indexInstability-1:indexInstability],freqsOfMode[mode][indexInstability-1:indexInstability],0)
    break
end

# Find divergence speed
divergenceSpeed = NaN
indicesNonOscillatoryInstability = [findfirst(x->x>0,problem[i].dampingsNonOscillatory) for i in eachindex(URange)]
indexDivergence = findfirst(!isnothing,indicesNonOscillatoryInstability)
divergenceSpeed = !isnothing(indexDivergence) ? URange[indexDivergence] : NaN
# Note: the value of the first dampingsNonOscillatory crosses zero at an airspeed between 37.2 and 37.3, but once it becomes positive, it disappears. So the divergence speed does match very closely that of the reference solution below

# Reference solution by Patil & Hodges & Cesnik: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)
flutterSpeedRef = 32.21
flutterFreqRef = 22.61
divergenceSpeedRef = 37.29

# Compute relative differences
ϵUf = flutterSpeed/flutterSpeedRef - 1
ϵFf = flutterFreq/flutterFreqRef - 1
ϵUd = divergenceSpeed/divergenceSpeedRef - 1

# Display current solution and comparison to reference
println("Flutter speed = $flutterSpeed m/s, flutter frequency = $flutterFreq rad/s, divergence speed = $divergenceSpeed m/s")
println("Respective relative differences: $ϵUf, $ϵFf, $ϵUd")

println("Finished SMWLinearFlutter.jl")