using AeroBeams, LinearAlgebra, LinearInterpolations

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Airfoil
airfoil = deepcopy(NACA0012)

# Altitude
h = 20e3

# Gravity
g = 9.80665

# Discretization
nElem = 16

# Stiffness factor
λ = 1

# Pitch angle
θ = 0*π/180

# Beam curvatures
k1 = 0.0
k2 = 0.0

# Model
SMWFlutter,_ = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,θ=θ*π/180,k1=k1,k2=k2,nElem=nElem,altitude=h,g=g,stiffnessFactor=λ)

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,displayStatus=false)

# Set number of vibration modes
nModes = 5

# Set airspeed range and initialize outputs
URange = collect(0:0.2:35)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(URange))
u1_of_x1 = Array{Vector{Float64}}(undef,length(URange))
u2_of_x1 = Array{Vector{Float64}}(undef,length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(URange))
α_of_x1 = Array{Vector{Float64}}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))
problem = Array{EigenProblem}(undef,length(URange))
u2_modeShapes = Array{Vector{Float64}}(undef,length(URange),nModes)
u3_modeShapes = Array{Vector{Float64}}(undef,length(URange),nModes)
p1_modeShapes = Array{Vector{Float64}}(undef,length(URange),nModes)

# Undeformed nodal and midpoint positions
x1_0 = vcat([vcat(SMWFlutter.beams[1].elements[e].r_n1[1],SMWFlutter.beams[1].elements[e].r_n2[1]) for e in 1:nElem]...)
x3_0 = vcat([vcat(SMWFlutter.beams[1].elements[e].r_n1[3],SMWFlutter.beams[1].elements[e].r_n2[3]) for e in 1:nElem]...)
x1_e = [SMWFlutter.beams[1].elements[e].x1 for e in 1:nElem]

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update velocity of basis A 
    set_motion_basis_A!(model=SMWFlutter,v_A=[0;U;0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=SMWFlutter,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],normalizeModeShapes=true)
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
    # Displacements over span
    u1_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
    u2_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[2],problem[i].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
    u3_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    u1e_of_x1 = [problem[i].elementalStatesOverσ[end][e].u[1] for e in 1:nElem]
    # Tip OOP displacement
    tip_u3[i] = problem[i].nodalStatesOverσ[end][nElem].u_n2[3]
    # Angle of attack over span
    α_of_x1[i] = [problem[i].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:nElem]
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1_of_x1[i]
    x3_def[i] = x3_0 .+ u3_of_x1[i]
    # Bending and torsional mode shapes
    for m in 1:nModes
        u2_modeShapes[i,m] = [problem[i].modeShapesAbs[m].elementalStates[e].u[2] for e in 1:nElem]
        u3_modeShapes[i,m] = [problem[i].modeShapesAbs[m].elementalStates[e].u[3] for e in 1:nElem]
        p1_modeShapes[i,m] = [problem[i].modeShapesAbs[m].elementalStates[e].p[3] for e in 1:nElem]
    end
end

# Frequencies and dampings after mode tracking
freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Update frequencies and dampings order on problem
for i in eachindex(URange)
    problem[i].frequenciesOscillatory = freqs[i]
    problem[i].dampingsOscillatory = damps[i]
end

# Update mode shapes' order
for i in eachindex(URange)
    u2_modeShapes[i,:] = u2_modeShapes[i,matchedModes[i]]
    u3_modeShapes[i,:] = u3_modeShapes[i,matchedModes[i]]
    p1_modeShapes[i,:] = p1_modeShapes[i,matchedModes[i]]
end

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

# Flutter speed and flutter frequency 
dampsOfMode = Array{Vector{Float64}}(undef,nModes)
freqsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOnsetSpeed = [Float64[] for _ in 1:nModes]
flutterOnsetFreq = [Float64[] for _ in 1:nModes]
flutterOnsetTipDisp = [Float64[] for _ in 1:nModes]
flutterOffsetSpeed = [Float64[] for _ in 1:nModes]
flutterOffsetFreq = [Float64[] for _ in 1:nModes]
flutterOffsetTipDisp = [Float64[] for _ in 1:nModes]
for mode in 1:nModes
    dampsOfMode[mode] = [damps[j][mode] for j in eachindex(URange)]
    freqsOfMode[mode] = [freqs[j][mode] for j in eachindex(URange)]
    # Flutter onset
    iOnset = 1 .+ findall(i -> dampsOfMode[mode][i] < 0 && dampsOfMode[mode][i+1] > 0, 1:length(dampsOfMode[mode])-1)
    if isempty(iOnset) || isempty(filter!(x->x!=1,iOnset))
        continue
    end
    for i in iOnset
        push!(flutterOnsetSpeed[mode],interpolate(dampsOfMode[mode][i-1:i],URange[i-1:i],0))
        push!(flutterOnsetFreq[mode],interpolate(dampsOfMode[mode][i-1:i],freqsOfMode[mode][i-1:i],0))
        push!(flutterOnsetTipDisp[mode],interpolate(dampsOfMode[mode][i-1:i],tip_u3[i-1:i],0))
    end
    # Flutter offset
    iOffset = 1 .+ findall(i -> dampsOfMode[mode][i] > 0 && dampsOfMode[mode][i+1] < 0, 1:length(dampsOfMode[mode])-1)
    if isempty(iOffset)
        continue
    end
    for i in iOffset
        push!(flutterOffsetSpeed[mode],interpolate(-dampsOfMode[mode][i-1:i],URange[i-1:i],0))
        push!(flutterOffsetFreq[mode],interpolate(-dampsOfMode[mode][i-1:i],freqsOfMode[mode][i-1:i],0))
        push!(flutterOffsetTipDisp[mode],interpolate(-dampsOfMode[mode][i-1:i],tip_u3[i-1:i],0))
    end
end
for mode in 1:nModes
    if isempty(flutterOnsetSpeed[mode])
        continue
    end
    println("Mode $mode: Flutter onset speed = $(flutterOnsetSpeed[mode]) m/s, flutter onset frequency = $(flutterOnsetFreq[mode]) rad/s")
end
for mode in 1:nModes
    if isempty(flutterOffsetSpeed[mode])
        continue
    end
    println("Mode $mode: Flutter offset speed = $(flutterOffsetSpeed[mode]) m/s, flutter offset frequency = $(flutterOffsetFreq[mode]) rad/s")
end

# Divergence speed
indicesNonOscillatoryInstability = [findfirst(x->x>0,problem[i].dampingsNonOscillatory) for i in eachindex(URange)]
indexDivergence = findfirst(!isnothing,indicesNonOscillatoryInstability)
divergenceSpeed = !isnothing(indexDivergence) ? URange[indexDivergence] : NaN
if isnan(divergenceSpeed) || divergenceSpeed == 0
    println("Divergence not found")
else
    println("Divergence speed = $divergenceSpeed m/s")
end

println("Finished SMWFlutter.jl")