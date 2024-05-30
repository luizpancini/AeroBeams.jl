using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Stiffness factor
λ = 1

# TF to include beam pods
beamPods = false

# Aerodynamic solver
aeroSolver = Indicial()

# Airspeed
U = 40*0.3048

# Option for mode tracking
modeTracking = true

# Number of modes
nModes = 10

# Set NR system solver for trim problem
relaxFactor = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor)

# Set payload range, and initialize outputs
PRange = collect(0:25:500)
untrackedFreqs = Array{Vector{Float64}}(undef,length(PRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(PRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(PRange))
freqs = Array{Vector{Float64}}(undef,length(PRange))
damps = Array{Vector{Float64}}(undef,length(PRange))

# Sweep payload
for (i,P) in enumerate(PRange)
    # Display progress
    println("Solving for payload = $P lb")
    # Model for trim problem
    heliosTrim,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
    # Set initial guess solution as previous known solution
    x0Trim = (i==1) ? zeros(0) : trimProblem.x
    # Create and solve trim problem
    global trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
    trimδ = trimProblem.x[end]
    # Model for eigen problem
    heliosEigen,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust)
    # Set initial solution as trim solution
    x0Eigen = trimProblem.x[1:end-2]
    # Create and solve eigen problem
    eigenProblem = create_EigenProblem(model=heliosEigen,x0=x0Eigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64])
    solve!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(eigenProblem.dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
end

# Apply mode tracking, if applicable
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(PRange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(PRange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(PRange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

# Plots
# ------------------------------------------------------------------------------
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(PRange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt12 = plot(xlabel="Payload [lb]", ylabel="Damping Ratio", ylims=[-1.0,0.25], legend=:bottomleft)
for mode in 1:nModes
    scatter!(PRange, modeDampingRatios[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/heliosFlutterPRange_1.pdf"))
# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-1.0,0.25])
for mode in 1:nModes
    scatter!(modeDampingRatios[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/heliosFlutterPRange_2.pdf"))

println("Finished heliosFlutterPRange.jl")