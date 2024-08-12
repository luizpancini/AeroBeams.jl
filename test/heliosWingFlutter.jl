using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Stiffness factor
λ = 1

# TF to include beam pod
beamPods = false

# Aerodynamic solver
aeroSolver = Indicial()

# Option for mode tracking
modeTracking = true

# Number of modes
nModes = 8

# Set system solver options
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set airspeed range, and initialize outputs
URange = collect(0.1:0.1:16.9)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model 
    _,_,heliosWing,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,airspeed=U)
    # plt = plot_undeformed_assembly(heliosWing)
    # display(plt)
    # Create and solve trim problem
    global eigenProblem = create_EigenProblem(model=heliosWing,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],systemSolver=NR)
    solve!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(eigenProblem.dampingsOscillatory,1e-8)
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
end

# Apply mode tracking, if applicable
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
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

# Plots
# ------------------------------------------------------------------------------
modeColors = get(colorschemes[:jet1], LinRange(0, 1, nModes))
lw = 2
ms = 3
relPath = "/test/outputs/figures/heliosWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Mode shapes
modesPlot = plot_mode_shapes(eigenProblem,scale=5,view=(30,30),frequencyLabel="frequency",save=true,savePath=string(relPath,"/heliosWingFlutter_modeShapes.pdf"))
display(modesPlot)
# V-g-f
gr()
plt11 = plot(ylabel="Frequency [rad/s]", ylims=[0,15])
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-2.0,0.25], legend=:bottomleft)
for mode in 1:nModes
    scatter!(URange, modeDampingRatios[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosWingFlutter_Vgf.pdf"))
# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-2.0,0.25] ,ylims=[0,15])
for mode in 1:nModes
    scatter!(modeDampingRatios[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosWingFlutter_rootlocus.pdf"))

println("Finished heliosWingFlutter.jl")