using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes

# Wing surface
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=flatPlate,c=chord,normSparPos=normSparPos)

# Set root angle
θ = π/180*1.3

# Wing beam
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIs = 0.75,0.1
ρIy = ρIs*EIy/EIz
ρIz = ρIs-ρIy
nElem = 16
∞ = 1e12
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
g = 9.807
h = 20e3
SMWFlutter = create_Model(name="SMWFlutter",beams=[wing],BCs=[clamp],gravityVector=[0;0;-g],altitude=h)

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set number of vibration modes
nModes = 5

# Set airspeed range and initialize outputs
URange = collect(15:0.25:33)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
u1_of_x1 = Array{Vector{Float64}}(undef,length(URange))
u2_of_x1 = Array{Vector{Float64}}(undef,length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))
problem = Array{EigenProblem}(undef,length(URange))
u2_modeShapes = Array{Vector{Float64}}(undef,length(URange),nModes)
u3_modeShapes = Array{Vector{Float64}}(undef,length(URange),nModes)
p1_modeShapes = Array{Vector{Float64}}(undef,length(URange),nModes)

# Undeformed nodal positions
x1 = vcat([vcat(SMWFlutter.beams[1].elements[e].x1_n1,SMWFlutter.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
x1_e = [SMWFlutter.beams[1].elements[e].x1 for e in 1:nElem]

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update velocity of basis A 
    set_motion_basis_A!(model=SMWFlutter,v_A=[0;U;0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=SMWFlutter,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64],normalizeModeShapes=true)
    solve!(problem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem[i].frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem[i].dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = problem[i].eigenvectorsOscillatoryCplx
    # Displacements over span
    u1_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
    u2_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[2],problem[i].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
    u3_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    # Deformed nodal positions
    x1_def[i] = x1 .+ u1_of_x1[i]
    x3_def[i] = u3_of_x1[i]
    # Bending and torsional mode shapes
    for m in 1:nModes
        u2_modeShapes[i,m] = [problem[i].modeShapesAbs[m].elementalStates[e].u[2] for e in 1:nElem]
        u3_modeShapes[i,m] = [problem[i].modeShapesAbs[m].elementalStates[e].u[3] for e in 1:nElem]
        p1_modeShapes[i,m] = [problem[i].modeShapesAbs[m].elementalStates[e].p[3] for e in 1:nElem]
    end
end

# Frequencies and dampings after mode tracking
freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Update mode shapes' order
for i in eachindex(URange)
    u2_modeShapes[i,:] = u2_modeShapes[i,matchedModes[i]]
    u3_modeShapes[i,:] = u3_modeShapes[i,matchedModes[i]]
    p1_modeShapes[i,:] = p1_modeShapes[i,matchedModes[i]]
end

# Separate frequencies and damping ratios by mode
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampingRatios[mode] = [damps[i][mode]/freqs[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end

# Flutter speed and flutter frequency 
dampsOfMode = Array{Vector{Float64}}(undef,nModes)
freqsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterSpeedOfMode = Array{Float64}(undef,nModes)
flutterFreqOfMode = Array{Float64}(undef,nModes)
for mode in 1:nModes
    dampsOfMode[mode] = [damps[j][mode] for j in eachindex(URange)]
    freqsOfMode[mode] = [freqs[j][mode] for j in eachindex(URange)]
    indexInstability = findfirst(x->x>0,dampsOfMode[mode])
    if isnothing(indexInstability)
        flutterSpeedOfMode[mode] = NaN
        flutterFreqOfMode[mode] = NaN
        continue
    end
    flutterSpeedOfMode[mode] = interpolate(dampsOfMode[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
    flutterFreqOfMode[mode] = interpolate(dampsOfMode[mode][indexInstability-1:indexInstability],freqsOfMode[mode][indexInstability-1:indexInstability],0)
end
if all(isnan,flutterSpeedOfMode)
    println("Flutter not found")
else
    flutterSpeed = minimum(filter(!isnan,flutterSpeedOfMode))
    flutterMode = findfirst(x->x==flutterSpeed,flutterSpeedOfMode)
    flutterFreq = flutterFreqOfMode[flutterMode]
    println("Flutter speed = $flutterSpeed m/s, flutter frequency = $flutterFreq rad/s (mode $flutterMode)")
end

# Divergence speed
indicesNonOscillatoryInstability = [findfirst(x->x>0,problem[i].dampingsNonOscillatory) for i in eachindex(URange)]
indexDivergence = findfirst(!isnothing,indicesNonOscillatoryInstability)
divergenceSpeed = !isnothing(indexDivergence) ? URange[indexDivergence] : NaN
if isnan(divergenceSpeed)
    println("Divergence not found")
else
    println("Divergence speed = $divergenceSpeed m/s")
end

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
for i in eachindex(URange)
    plot!(x1_def[i]/L, x3_def[i]/L, c=colors[i], lw = lw, label=false)
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutter_1.pdf"))
# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-0.5,0.1])
for mode in 1:nModes
    plot!(modeDampingRatios[mode], modeFrequencies[mode], c=modeColors[mode], lw=lw, label="Mode $mode")
end
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutter_2.pdf"))
# V-g-f
plt31 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], lw=lw,  label=false)
end
plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.05])
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutter_3.pdf"))
# Mode shapes at specified velocity
U2plot = URange[1]
ind = findfirst(x->x==U2plot,URange)
plt4 = plot(xlabel="\$x_1/L\$", ylabel="\$u_2\$", title="Chordwise bending mode shape at U=$(URange[ind]) m/s")
plt5 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3\$", title="Flapwise bending mode shape at U=$(URange[ind]) m/s")
plt6 = plot(xlabel="\$x_1/L\$", ylabel="\$p_1\$", title="Torsional mode shape at U=$(URange[ind]) m/s")
for m in 1:nModes
    plot!(plt4,x1_e/L, u2_modeShapes[ind,m], c=modeColors[m], lw = lw, label="Mode $m")
    plot!(plt5,x1_e/L, u3_modeShapes[ind,m], c=modeColors[m], lw = lw, label="Mode $m")
    plot!(plt6,x1_e/L, p1_modeShapes[ind,m], c=modeColors[m], lw = lw, label="Mode $m")
end
display(plt4)
display(plt5)
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutter_4.pdf"))

println("Finished SMWFlutter.jl")