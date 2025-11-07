using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Flight altitude
h = 20e3*0.3048

# Airspeed range
URange = collect(60:1:180)

# Number of vibration modes
nModes = 5

# Flag to update airfoil parameters (lift-slope)
updateAirfoilParameters = true

# System solver
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter)

# Pre-allocate memory and initialize output arrays
trimProblem = Array{TrimProblem}(undef,length(URange))
eigenProblem = Array{EigenProblem}(undef,length(URange))
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed range
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    BWBtrim = first(create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true,updateAirfoilParameters=updateAirfoilParameters))
    # Initial guess for trim problem
    x0Trim = i == 1 ? zeros(0) : trimProblem[i-1].x
    # Solve the trim problem
    trimProblem[i] = create_TrimProblem(model=BWBtrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem[i])
    # Trim outputs
    trimAoA[i] = trimProblem[i].aeroVariablesOverσ[end][BWBtrim.beams[3].elementRange[1]].flowAnglesAndRates.αₑ
    trimThrust[i] = trimProblem[i].x[end-1]*BWBtrim.forceScaling 
    trimδ[i] = trimProblem[i].x[end]
    println("Trim AoA = $(round(trimAoA[i]*180/π,digits=2)) deg, trim thrust = $(round(trimThrust[i],digits=2)) N, trim δ = $(round(trimδ[i]*180/π,digits=2)) deg")
    # Wing model for eigenproblem
    _,rightWingModel = create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElev=trimδ[i],wingRootAoA=trimAoA[i],thrust=trimThrust[i],updateAirfoilParameters=updateAirfoilParameters)
    # Initial guess for eigenproblem
    x0Eig = i == 1 ? zeros(0) : eigenProblem[i-1].x
    # Solve eigenproblem
    eigenProblem[i] = create_EigenProblem(model=rightWingModel,nModes=nModes,frequencyFilterLimits=[2e0,Inf],systemSolver=NR,x0=x0Eig)
    solve!(eigenProblem[i])
    # Extract eigen-solution
    untrackedFreqs[i] = eigenProblem[i].frequenciesOscillatory
    untrackedDamps[i] = eigenProblem[i].dampingsOscillatory
    untrackedEigenvectors[i] = eigenProblem[i].eigenvectorsOscillatoryCplx
end

# Mode tracking
freqs,damps,_,matchedModes = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end

# Flutter onset speeds and frequencies
flutterOnsetSpeedOfMode = fill(NaN, nModes)
flutterOnsetFreqOfMode = fill(NaN, nModes)
for mode in 1:nModes
    iOnset = findfirst(j -> modeDampings[mode][j] < 0 && modeDampings[mode][j+1] > 0, 1:length(URange)-1)
    flutterOnsetSpeedOfMode[mode] = isnothing(iOnset) ? Inf : interpolate(modeDampings[mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
    flutterOnsetFreqOfMode[mode] = isnothing(iOnset) ? Inf : interpolate(modeDampings[mode][iOnset:iOnset+1],modeFrequencies[mode][iOnset:iOnset+1],0)
end
flutterModes = findall(!isinf, flutterOnsetSpeedOfMode)
flutterOnsetSpeedsAll = filter(!isinf,flutterOnsetSpeedOfMode)
flutterOnsetFreqsAll = filter(!isinf,flutterOnsetFreqOfMode)

flutterMode = flutterModes[argmin(flutterOnsetSpeedsAll)]
flutterOnsetSpeed = minimum(flutterOnsetSpeedsAll)
flutterSpeedInd = findfirst(x-> x>flutterOnsetSpeed, URange)
for (mode,speed,freq) in zip(flutterModes,flutterOnsetSpeedsAll,flutterOnsetFreqsAll)
    println("Mode $mode: flutter speed = $(round(speed,sigdigits=4)) m/s, flutter frequency = $(round(freq/(2π),digits=2)) Hz")
end

println("Finished BWBmatchedWingFlutter.jl")