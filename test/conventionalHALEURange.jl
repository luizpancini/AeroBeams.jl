using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Stiffness factor for the aircraft's wing
λ = 1e0

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10

# System solver for trim problem
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=false)

# Set number of vibration modes
nModes = 15

# Set airspeed range and initialize outputs
URange = collect(20:0.5:35)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementID=1,localNode=1,ku=ku,kp=kp)
spring2 = create_Spring(elementID=nElemTailBoom,localNode=2,ku=ku,kp=kp)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    conventionalHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true)
    # Add springs
    add_springs_to_beam!(tailBoom,springs=[spring1,spring2])
    # Update model
    conventionalHALEtrim.skipValidationMotionBasisA = true
    update_model!(conventionalHALEtrim)
    # Set initial guess solution as previous known solution
    x0Trim = i == 1 ? zeros(0) : trimProblem.x
    # Create and trim problem
    global trimProblem = create_TrimProblem(model=conventionalHALEtrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimAoA[i] = trimProblem.flowVariablesOverσ[end][div(nElemWing,2)].αₑ
    trimThrust[i] = stabilizersAero ? trimProblem.x[end-1]*trimProblem.model.forceScaling : trimProblem.x[end]*trimProblem.model.forceScaling
    trimδ[i] = stabilizersAero ? trimProblem.x[end] : 0
    println("Trim AoA = $(trimAoA[i]*180/π), trim thrust = $(trimThrust[i]), trim δ = $(trimδ[i]*180/π)")
    # Model for eigen problem
    conventionalHALEeigen,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ[i],thrust=trimThrust[i])
    # Add springs
    add_springs_to_beam!(tailBoom,springs=[spring1,spring2])
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
# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[20,35], ylims=[0,15])
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEURange_AoA.pdf"))
# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,35])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEURange_thrust.pdf"))
# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,35])
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEURange_delta.pdf"))
# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-16.0,1.0], ylims=[0,50])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEURange_rootlocus.pdf"))
# Root locus (zoom)
plt5 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-1,1], ylims=[0,10])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEURange_rootlocuszoom.pdf"))


println("Finished conventionalHALEURange.jl")