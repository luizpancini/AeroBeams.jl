using AeroBeams

# Airspeed range [m/s]
URange = collect(20:0.5:35)

# Stiffness factor for the aircraft's wing
λ = 1e0

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Discretization
nElemWing = 40
nElemTailBoom = 5
nElemHorzStabilizer = 4

# Attachment spring stiffness
μ = 1e-1

# System solver for trim problem
NR = create_NewtonRaphson(trimSolverMethod=:LM,displayStatus=false)

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :arnoldi
nev = 140

# Initialize outputs
eigenProblem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

# Initial guess for trim problem
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    cHALEtrim,_ = create_cHALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,withAttachmentSprings=true,μ=μ)
    # Create and trim problem
    trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR,x0=x0)
    solve!(trimProblem) 
    # Extract trim variables
    trimSol = trimProblem.trimSolution
    aeroVars = trimProblem.loadstepHistory.aeroVariables
    trimAoA[i] = aeroVars[end][div(nElemWing,2)].unsteadyKin.αₑ
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = stabilizersAero ? trimSol.trimδ[1] : NaN
    println("Trim AoA = $(round(rad2deg(trimAoA[i]), digits=2)) deg, trim thrust = $(round(trimThrust[i], digits=2)) N, trim δ = $(round(rad2deg(trimδ[i]), digits=2)) deg")
    # Model for eigen problem
    cHALEeigen,_ = create_cHALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ[i],thrust=trimThrust[i],withAttachmentSprings=true,μ=μ)
    # Create and solve eigen problem
    eigenProblem[i] = create_EigenProblem(model=cHALEeigen,eigenmethod=eigenmethod,nev=nev,refTrimProblem=trimProblem)
    solve_eigen!(eigenProblem[i])
    # Eigensolution
    eigensols[i] = eigenProblem[i].eigensolution
    # Set initial guess solution for next iteration
    global x0 = trimProblem.systemArrays.x
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, fmax=20, maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALEURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at last airspeed
modesPlot = plot_mode_shapes(eigenProblem[end],modeShapes=modeShapesOsc[:,end],modeNumbers=collect(1:nModes),scale=10,view=(30,30),legendPos=(0.3,0.98),save=true,savePath=string(relPath,"/cHALEURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = cgrad(:rainbow, length(URange), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[20,35], ylims=[0,15])
plot!(URange, rad2deg.(trimAoA), c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALEURange_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,35])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALEURange_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,35])
plot!(URange, rad2deg.(trimδ), c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/cHALEURange_delta.pdf"))

# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-16.0,1.0], ylims=[0,12])
for mode in 1:nModes
    scatter!(dampsOsc[mode,:]./freqsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/cHALEURange_rootlocus.pdf"))

# Root locus (zoom)
plt5 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-4,0.1], ylims=[0,12])
for mode in 1:nModes
    scatter!(dampsOsc[mode,:]./freqsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt5)
savefig(string(absPath,"/cHALEURange_rootlocuszoom.pdf"))

println("Finished cHALEURange.jl")