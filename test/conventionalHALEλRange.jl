using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Mode tracking option
modeTracking = false

# Airspeed
U = 25

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=100)

# Set number of vibration modes
nModes = 12

# Set stiffness factor range and initialize outputs
λRange = [50,40,30,20,18,16,14,12,10,9,8,7,6,5,4,3,2.5,2,1.75,1.5,1.25,1.1,1]
untrackedFreqs = Array{Vector{Float64}}(undef,length(λRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(λRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(λRange))
freqs = Array{Vector{Float64}}(undef,length(λRange))
damps = Array{Vector{Float64}}(undef,length(λRange))

# Sweep airspeed
for (i,λ) in enumerate(λRange)
    println("Solving for λ = $λ")
    # Model
    conventionalHALE,_ = create_conventional_HALE(aeroSolver=Inflow(),nElemWing=16,nElemTailBoom=10,nElemHorzStabilizer=8,stiffnessFactor=λ,∞=1e12,stabilizersAero=true,includeVS=false,wingCd0=1e-2)
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
    freqs,damps,_,matchedModes = mode_tracking(λRange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(λRange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(λRange)]
end

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 4
xtick_positions = [1,2,5,10,20,50]
xtick_labels = ["1","2","5","10","20","50"]
ytick_positions = [0,10,20,50]
ytick_labels = ["0","10","20","50"]
# Frequency and damping evolution
plt31 = plot(ylabel="Frequency [rad/s]", xlims=[0.9,maximum(λRange)], ylims=[0,50], xticks=(xtick_positions, xtick_labels), yticks=(ytick_positions, ytick_labels), xscale=:log10, yscale=:log10)
for mode in 1:nModes
    scatter!(λRange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0,  label=false)
end
plt32 = plot(xlabel="Stiffness factor", ylabel="Damping [1/s]", xlims=[0.9,maximum(λRange)], ylims=[-5,2], xscale=:log10, xticks=(xtick_positions, xtick_labels))
for mode in 1:nModes
    scatter!(λRange, modeDampings[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEλRange.pdf"))

println("Finished conventionalHALEλRange.jl")