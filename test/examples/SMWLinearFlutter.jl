using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes

# Wing surface
airfoil = deepcopy(flatPlate)
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos)

# Wing beam
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIs = 0.75,0.1
ρIy = ρIs*EIy/EIz
ρIz = ρIs-ρIy
nElem = 16
∞ = 1e12
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)],aeroSurface=surf)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
g = 9.80665
h = 20e3
SMWLinearFlutter = create_Model(name="SMWLinearFlutter",beams=[wing],BCs=[clamp],gravityVector=[0;0;-g],altitude=h)

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

# Plots
# ------------------------------------------------------------------------------
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
relPath = "/test/outputs/figures/SMWLinearFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Plot mode shapes
modesPlot = plot_mode_shapes(problem[end],scale=5,view=(30,30),legendPos=:best,frequencyLabel="frequency",save=true,savePath=string(relPath,"/SMWLinearFlutter_modeShapes.pdf"))
display(modesPlot)
# V-g-f
gr()
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], mc=modeColors[mode], ms = ms, msw=0, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.2,0.1])
for mode in 1:nModes
    scatter!(URange, modeDampingRatios[mode], mc=modeColors[mode], ms = ms, msw=0, label=false)
    scatter!([NaN], [NaN], mc=modeColors[mode], ms = ms, msw=0, label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWLinearFlutter_Vgf.pdf"))

println("Finished SMWLinearFlutter.jl")