using AeroBeams

# Stiffness range
λRange = [1,2,20]

# Airspeed range [m/s]
URange = collect(20:0.5:35)

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Discretization
nElemWing = 40
nElemTailBoom = 1
nElemHorzStabilizer = 2

# Attachment spring stiffness
μ = 1e-1

# Aerodynamic solver
aeroSolver = Indicial()

# System solver for trim problem
NR = create_NLsolver(trimSolverMethod=:LQ)

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :arnoldi
nev = 150

# Initialize outputs
eigenProblem = Array{EigenProblem}(undef,length(λRange),length(URange))
eigensols = Array{EigenSolution}(undef,length(λRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(λRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(λRange))
modeShapesOsc = Array{Matrix{Union{Nothing,ModeShape}}}(undef,length(λRange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving trim problem for λ = $λ, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_ = create_cHALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e10,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,withAttachmentSprings=true,μ=μ)
        # Create and trim problem
        trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR,x0=x0)
        solve!(trimProblem)
        # Extract trim variables
        trimSol = trimProblem.trimSolution
        aeroVars = trimProblem.loadstepHistory.aeroVariables
        trimAoA = aeroVars[end][div(nElemWing,2)].unsteadyKin.αₑ
        trimThrust = trimSol.trimLoads[1]
        trimδ = stabilizersAero ? trimSol.trimδ[1] : NaN
        println("Trim AoA = $(round(rad2deg(trimAoA), digits=2)) deg, trim thrust = $(round(trimThrust, digits=2)) N, trim δ = $(round(rad2deg(trimδ), digits=2)) deg")
        # Model for eigen problem
        cHALEeigen,_ = create_cHALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ,thrust=trimThrust,withAttachmentSprings=true,μ=μ)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=cHALEeigen,eigenmethod=eigenmethod,nev=nev,refTrimProblem=trimProblem)
        solve_eigen!(eigenProblem[i,j])
        # Eigensolution
        eigensols[i,j] = eigenProblem[i,j].eigensolution
        # Set initial guess solution for next iteration
        x0 = trimProblem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,modeShapes,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],modeShapesOsc[i],_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, fmax=20, maxTracks=nModes)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALELURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes of flexible aircraft at lowest airspeed
ind = 1
modesPlot = plot_mode_shapes(eigenProblem[1,ind],modeShapes=modeShapesOsc[1][:,ind],modeNumbers=collect(1:nModes),scale=5,view=(30,30),legendPos=(0.3,0.99),save=true,savePath=string(relPath,"/cHALELURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = cgrad(:rainbow, length(λRange), categorical=true)
lw = 2
ms = 3
msw = 0
lstyle = [:solid, :dash, :dot]
mshape = [:circle, :star, :utriangle]
gr()

# Root locus
plt0 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", ylims=[0,20])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=msw, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in axes(freqsOsc[i], 1)
        scatter!(dampsOsc[i][mode,:]./freqsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
display(plt0)
savefig(string(absPath,"/cHALELURange_rootlocus.pdf"))

# Root locus
plt1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,0], ylims=[0,10])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=msw, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in axes(freqsOsc[i], 1)
        scatter!(dampsOsc[i][mode,:]./freqsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/cHALELURange_rootlocus2.pdf"))

# Root locus (zoom)
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.3,0.05], ylims=[0,0.5])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=msw, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in axes(freqsOsc[i], 1)
        scatter!(dampsOsc[i][mode,:]./freqsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/cHALELURange_rootlocuszoom.pdf"))

println("Finished cHALELURange.jl")