using AeroBeams

# Stiffness factor range
λRange = [50,40,30,25,20,18,16,14,12,10,9,8,7,6,5,4,3,2.5,2,1.75,1.5,1.25,1.1,1]

# Airspeed [m/s]
U = 25

# Aerodynamic solver
aeroSolver = Indicial()

# Discretization
nElemWing = 40
nElemTailBoom = 5
nElemHorzStabilizer = 4

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Attachment spring stiffness
μ = 1e-1

# System solver for trim problem
NR = create_NewtonRaphson(trimSolverMethod=:LQ,displayStatus=false)

# Number of modes to track
nModes = 20

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 120

# Initialize outputs
eigenProblem = Vector{EigenProblem}(undef,length(λRange))
trimAoA = Array{Float64}(undef,length(λRange))
trimThrust = Array{Float64}(undef,length(λRange))
trimδ = Array{Float64}(undef,length(λRange))
eigensols = Array{EigenSolution}(undef,length(λRange))

# Initial guess for trim problem
x0 = zeros(0)

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    println("Solving for λ = $λ")
    # Model for trim problem
    cHALEtrim,_ = create_cHALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,withAttachmentSprings=true,μ=μ)
    # Create and trim problem
    trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR,x0=x0)
    solve!(trimProblem)
    # Extract trim variables
    trimSol = trimProblem.trimSolution
    aeroVars = trimProblem.loadstepHistory.aeroVariables
    trimAoA[i] = rad2deg(aeroVars[end][div(nElemWing,2)].unsteadyKin.αₑ)
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = stabilizersAero ? trimSol.trimδ[1] : NaN
    println("Trim AoA = $(round(trimAoA[i], digits=2)) deg, trim thrust = $(round(trimThrust[i], digits=2)) N, trim δ = $(round(rad2deg(trimδ[i]), digits=2)) deg")
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
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALELambdaRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes of nominal-stiffness aircraft 
modesPlot = plot_mode_shapes(eigenProblem[end],modeShapes=modeShapesOsc[:,end],modeNumbers=collect(1:nModes),scale=10,view=(30,30),legendPos=(-0.3,0.99),save=true,savePath=string(relPath,"/cHALELambdaRange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = cgrad(:rainbow, length(λRange), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 4
xtick_positions = [1,2,5,10,20,50]
xtick_labels = ["1","2","5","10","20","50"]
ytick_positions = [0,10,20,50]
ytick_labels = ["0","10","20","50"]
gr()

# Trim root angle of attack vs stiffness factor
plt1 = plot(xlabel="Stiffness factor", ylabel="Trim root AoA [deg]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimAoA, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALELambdaRange_AoA.pdf"))

# Trim propeller force vs stiffness factor
plt2 = plot(xlabel="Stiffness factor", ylabel="Trim thrust [N]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALELambdaRange_thrust.pdf"))

# Trim elevator deflection vs stiffness factor
plt3 = plot(xlabel="Stiffness factor", ylabel="Trim elevator deflection [deg]", xlims=[minimum(λRange),maximum(λRange)], xscale=:log10, xticks=(xtick_positions, xtick_labels))
plot!(λRange, rad2deg.(trimδ), c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/cHALELambdaRange_delta.pdf"))

# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,1], ylims=[0,300], yscale=:log10)
for mode in 1:nModes
    scatter!(dampsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/cHALELambdaRange_rootlocus.pdf"))

# Frequency and damping evolution
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[0.9,maximum(λRange)], ylims=[0,50], xticks=(xtick_positions, xtick_labels), yticks=(ytick_positions, ytick_labels), xscale=:log10, yscale=:log10)
for mode in 1:nModes
    plot!(λRange, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt52 = plot(xlabel="Stiffness factor", ylabel="Damping [1/s]", xlims=[0.9,maximum(λRange)], ylims=[-5,2], xscale=:log10, xticks=(xtick_positions, xtick_labels))
for mode in 1:nModes
    plot!(λRange, dampsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/cHALELambdaRange_lambda-g-f.pdf"))

println("Finished cHALELambdaRange.jl")