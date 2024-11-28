using AeroBeams, LinearInterpolations

# Goland wing parameters taken from PALACIOS & EPUREANU - An Intrinsic Description of the Nonlinear Aeroelasticity of Very Flexible Wings - (2011) (https://doi.org/10.2514/6.2011-1917)

# Aerodynamic solver
aeroSolver = Inflow(6)

# Derivation method
derivationMethod = AD()

# Gravity
g = 0

# Altitude
altitude = 1.867e3
atmosphere = standard_atmosphere(altitude)

# Pitch angle
θ = 0*π/180

# Number of modes
nModes = 2

# Airspeed range
URange = collect(20:5:160)

# Tip loss factor
τ = Inf

# Aerodynamic surface
airfoil = deepcopy(flatPlate)
chord = 1.8288
normSparPos = 0.33
normCGPos = 0.43
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=true,tipLossDecayFactor=τ,updateAirfoilParameters=true)

# Wing beam
L = 6.096
EIy = 9.77e6
GJ = 0.99e6
ρA = 35.71
ρIs = 8.64
e2 = -(normCGPos-normSparPos)*chord
∞ = 1e12
nElem = 30
beam = create_Beam(name="wingBeam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy)],I=[inertia_matrix(ρA=ρA,ρIs=ρIs,e2=e2)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
GolandWing = create_Model(name="GolandWing",beams=[beam],BCs=[clamp],gravityVector=[0;0;g],atmosphere=atmosphere)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
dampingsNonOscillatory = Array{Vector{Float64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
problem = Array{EigenProblem}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    ## Update airspeed on model
    set_motion_basis_A!(model=GolandWing,v_A=[0;U;0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=GolandWing,nModes=nModes)
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
    dampingsNonOscillatory[i] = problem[i].dampingsNonOscillatory
end

# Apply mode tracking
freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
end

## Flutter speed and flutter frequency 
flutterOnsetSpeed = [Float64[] for _ in 1:nModes]
flutterOnsetFreq = [Float64[] for _ in 1:nModes]
flutterOffsetSpeed = [Float64[] for _ in 1:nModes]
flutterOffsetFreq = [Float64[] for _ in 1:nModes]
for mode in 1:nModes
    ## Flutter onset
    iOnset = 1 .+ findall(i -> modeDampings[mode][i] < 0 && modeDampings[mode][i+1] > 0, 1:length(modeDampings[mode])-1)
    if isempty(iOnset) || isempty(filter!(x->x!=1,iOnset))
        continue
    end
    for i in iOnset
        push!(flutterOnsetSpeed[mode],interpolate(modeDampings[mode][i-1:i],URange[i-1:i],0))
        push!(flutterOnsetFreq[mode],interpolate(modeDampings[mode][i-1:i],modeFrequencies[mode][i-1:i],0))
    end
    ## Flutter offset
    iOffset = 1 .+ findall(i -> modeDampings[mode][i] > 0 && modeDampings[mode][i+1] < 0, 1:length(modeDampings[mode])-1)
    if isempty(iOffset)
        continue
    end
    for i in iOffset
        push!(flutterOffsetSpeed[mode],interpolate(-modeDampings[mode][i-1:i],URange[i-1:i],0))
        push!(flutterOffsetFreq[mode],interpolate(-modeDampings[mode][i-1:i],modeFrequencies[mode][i-1:i],0))
    end
end

# Show flutter speed and frequency
flutterSpeedAll = vcat(flutterOnsetSpeed...)
flutterFreqAll = vcat(flutterOnsetFreq...)
ind = argmin(flutterSpeedAll)
flutterSpeed = flutterSpeedAll[ind]
flutterFreq = flutterFreqAll[ind]
println("Flutter speed = $(flutterSpeed) m/s")
println("Flutter frequency = $(flutterFreq) rad/s")

println("Finished GolandWingFlutter.jl")