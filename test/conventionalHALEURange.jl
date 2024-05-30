using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Mode tracking option
modeTracking = true

# Stiffness factor for the aircraft's wing
λ = 5.0

# Model
conventionalHALE,_ = create_conventional_HALE(aeroSolver=Indicial(),nElemWing=20,nElemTailBoom=10,nElemHorzStabilizer=10,stiffnessFactor=λ,∞=1e12,stabilizersAero=true,includeVS=false,wingCd0=1e-2)

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=100)

# Set number of vibration modes
nModes = 8

# Set airspeed range and initialize outputs
URange = collect(15:0.5:35)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update velocity of basis A 
    set_motion_basis_A!(model=conventionalHALE,v_A=[0;U;0])
    # Create and solve problem
    problem = create_EigenProblem(model=conventionalHALE,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
    solve!(problem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem.dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = problem.eigenvectorsOscillatoryCplx
end

# Frequencies and dampings after mode tracking
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
# V-g-f
plt31 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,20*λ])
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], lw=lw,  label=false)
end
plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=[URange[1],URange[end]], ylims=[-5,1])
for mode in 1:nModes
    plot!(URange, modeDampings[mode], c=modeColors[mode], lw=lw, label="Mode $mode")
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEURange_3.pdf"))

println("Finished conventionalHALEURange.jl")