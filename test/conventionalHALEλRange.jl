using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Mode tracking option
modeTracking = false

# Airspeed
U = 25

# Aerodynamic solver
aeroSolver = Indicial()

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# System solver for trim problem
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=false)

# Set number of vibration modes
nModes = 12

# Set stiffness factor range and initialize outputs
λRange = [50,40,30,25,20,18,16,14,12,10,9,8,7,6,5,4,3,2.5,2,1.75,1.5,1.25,1.1,1]
trimAoA = Array{Float64}(undef,length(λRange))
trimThrust = Array{Float64}(undef,length(λRange))
trimδ = Array{Float64}(undef,length(λRange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(λRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(λRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(λRange))
freqs = Array{Vector{Float64}}(undef,length(λRange))
damps = Array{Vector{Float64}}(undef,length(λRange))

# Attachment springs
μ = 1e-1
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=ku,kp=kp)

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    println("Solving for λ = $λ")
    # Model for trim problem
    conventionalHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true)
    # Add springs
    add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
    # Update model
    conventionalHALEtrim.skipValidationMotionBasisA = true
    update_model!(conventionalHALEtrim)
    # Set initial guess solution as previous known solution
    x0Trim = i == 1 ? zeros(0) : trimProblem.x
    # Create and trim problem
    global trimProblem = create_TrimProblem(model=conventionalHALEtrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimAoA[i] = trimProblem.aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
    trimThrust[i] = stabilizersAero ? trimProblem.x[end-1]*trimProblem.model.forceScaling : trimProblem.x[end]*trimProblem.model.forceScaling
    trimδ[i] = stabilizersAero ? trimProblem.x[end] : 0
    println("Trim AoA = $(trimAoA[i]*180/π), trim thrust = $(trimThrust[i]), trim δ = $(trimδ[i]*180/π)")
    # Model for eigen problem
    conventionalHALEeigen,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ[i],thrust=trimThrust[i])
    # Add springs
    add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
    # Update model
    conventionalHALEeigen.skipValidationMotionBasisA = true
    update_model!(conventionalHALEeigen)
    # Create and solve eigen problem
    eigenProblem = create_EigenProblem(model=conventionalHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],jacobian=trimProblem.jacobian[1:end,1:end-trimProblem.model.nTrimVariables],inertia=trimProblem.inertia)
    solve_eigen!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(eigenProblem.dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
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
# Trim root angle of attack vs stiffness factor
plt1 = plot(xlabel="Stiffness factor", ylabel="Trim root AoA [deg]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEλRange_AoA.pdf"))
# Trim propeller force vs stiffness factor
plt2 = plot(xlabel="Stiffness factor", ylabel="Trim thrust [N]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEλRange_thrust.pdf"))
# Trim elevator deflection vs stiffness factor
plt3 = plot(xlabel="Stiffness factor", ylabel="Trim elevator deflection [deg]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEλRange_delta.pdf"))
# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,1], ylims=[0,300], yscale=:log10)
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEλRange_rootlocus.pdf"))
# Frequency and damping evolution
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[0.9,maximum(λRange)], ylims=[0,50], xticks=(xtick_positions, xtick_labels), yticks=(ytick_positions, ytick_labels), xscale=:log10, yscale=:log10)
for mode in 1:nModes
    scatter!(λRange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0,  label=false)
end
plt52 = plot(xlabel="Stiffness factor", ylabel="Damping [1/s]", xlims=[0.9,maximum(λRange)], ylims=[-5,2], xscale=:log10, xticks=(xtick_positions, xtick_labels))
for mode in 1:nModes
    scatter!(λRange, modeDampings[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEλRange_lambda-g-f.pdf"))

println("Finished conventionalHALEλRange.jl")