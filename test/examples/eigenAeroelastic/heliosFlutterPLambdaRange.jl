using AeroBeams

# Stiffness factor range
λRange = [1,2,100]

# Payload range [lb]
PRange = collect(0:10:500)

# Airspeed [m/s]
U = feet2meters(40)

# Wing airfoil
wingAirfoil = deepcopy(HeliosWingAirfoil)

# Option for reduced chord
reducedChord = false

# Option to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Attachment spring stiffness
μ = 1e-1

# Trim system solver
NR = create_NewtonRaphson(trimSolverMethod=:LQ,displayStatus=false)

# Number of modes to track
nModes = 12

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :arnoldi
nev = 200

# Initialize outputs
eigenProblem = Array{EigenProblem}(undef,length(λRange),length(PRange))
eigensols = Array{EigenSolution}(undef,length(λRange),length(PRange))
freqsOsc = Array{Matrix{Float64}}(undef,length(λRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(λRange))
modeShapesOsc = Array{Matrix{Union{Nothing,ModeShape}}}(undef,length(λRange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Solving for λ = $λ, payload = $P lb")
        # Model for trim problem
        heliosTrim,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,payloadOnWing=payloadOnWing,withAttachmentSprings=true,μ=μ)
        # Create and solve trim problem
        trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR,x0=x0)
        solve!(trimProblem)
        # Trim outputs
        trimSol = trimProblem.trimSolution
        aeroVars = trimProblem.loadstepHistory.aeroVariables
        trimAoA = aeroVars[end][midSpanElem].unsteadyKin.αₑ
        trimThrust = trimSol.trimLoads[1]
        trimδ = trimSol.trimδ[1]
        println("Trim AoA = $(round(rad2deg(trimAoA), digits=2)) deg, trim thrust = $(round(trimThrust, digits=2)) N, trim δ = $(round(rad2deg(trimδ), digits=2)) deg")
        # Model for eigen problem
        heliosEigen,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,reducedChord=reducedChord,payloadOnWing=payloadOnWing,withAttachmentSprings=true,μ=μ)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=heliosEigen,eigenmethod=eigenmethod,nev=nev,refTrimProblem=trimProblem)
        solve_eigen!(eigenProblem[i,j])
        # Eigensolution
        eigensols[i,j] = eigenProblem[i,j].eigensolution
        # Set initial guess solution for next iteration
        x0 = trimProblem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,modeShapes,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],modeShapesOsc[i],_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, maxTracks=nModes)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/heliosFlutterPLambdaRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes of flexible aircraft at highest payload
ind = length(PRange)
modesPlot = plot_mode_shapes(eigenProblem[1,ind],element2centralize=15,scale=10,view=(30,30),legendPos=:outerright,modeShapes=modeShapesOsc[1][:,ind],modeNumbers=collect(1:nModes),save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = cgrad(:rainbow, length(λRange), categorical=true)
lw = 2
ms = 3
lstyle = [:solid, :dash, :dot]
mshape = [:circle, :star, :utriangle]
gr()

# Root locus 
plt1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(dampsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/heliosFlutterPLambdaRange_rootlocus.pdf"))

# Root locus (phugoid zoom)
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.15,0.15], ylims=[0,0.8])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(dampsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/heliosFlutterPLambdaRange_Phzoom.pdf"))

# Root locus (dutch roll zoom)
plt3 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.6,0], ylims=[0,0.6])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(dampsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/heliosFlutterPLambdaRange_DRzoom.pdf"))

# Root locus (short-period zoom)
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-4.5,-2.5], ylims=[0,5])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(dampsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt4)
savefig(string(absPath,"/heliosFlutterPLambdaRange_SPzoom.pdf"))

println("Finished heliosFlutterPLambdaRange.jl")