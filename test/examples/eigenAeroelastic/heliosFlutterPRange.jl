using AeroBeams

# Payload range [lb]
PRange = collect(0:20:500)

# Airspeed [m/s]
U = feet2meters(40)

# Stiffness factor
λ = 1

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

# System solver for trim problem
NR = create_NewtonRaphson(trimSolverMethod=:LQ,displayStatus=false)

# Number of modes to track
nModes = 12

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = aeroSolver isa DynamicStallSolver ? 250 : 150

# Initialize outputs
eigenProblem = Array{EigenProblem}(undef,length(PRange))
eigensols = Array{EigenSolution}(undef,length(PRange))
trimAoA = Array{Float64}(undef,length(PRange))
trimThrust = Array{Float64}(undef,length(PRange))
trimδ = Array{Float64}(undef,length(PRange))

# Initial guess solution for trim problem
x0 = zeros(0)

# Sweep payload
for (i,P) in enumerate(PRange)
    # Display progress
    println("Solving for payload = $P lb")
    # Model for trim problem
    heliosTrim,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,payloadOnWing=payloadOnWing,withAttachmentSprings=true,μ=μ)
    # Create and solve trim problem
    trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR,x0=x0)
    solve!(trimProblem)
    # Extract trim variables
    trimSol = trimProblem.trimSolution
    aeroVars = trimProblem.loadstepHistory.aeroVariables
    trimAoA[i] = aeroVars[end][midSpanElem].unsteadyKin.αₑ
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = trimSol.trimδ[1]
    println("Trim AoA = $(round(rad2deg(trimAoA[i]), digits=2)) deg, trim thrust = $(round(trimThrust[i], digits=2)) N, trim δ = $(round(rad2deg(trimδ[i]), digits=2)) deg")
    # Model for eigen problem
    heliosEigen,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ[i],thrust=trimThrust[i],reducedChord=reducedChord,payloadOnWing=payloadOnWing,withAttachmentSprings=true,μ=μ)
    # Create and solve eigen problem
    eigenProblem[i] = create_EigenProblem(model=heliosEigen,eigenmethod=eigenmethod,nev=nev,refTrimProblem=trimProblem)
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
relPath = "/test/outputs/figures/heliosFlutterPRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at highest payload
modesPlot = plot_mode_shapes(eigenProblem[end],backendSymbol=:gr,element2centralize=15,scale=10,modeShapes=modeShapesOsc[:,end],modeNumbers=collect(1:nModes),view=(30,30),legendPos=:outerright,save=true,savePath=string(relPath,"/heliosFlutterPRange_modeShapes.pdf"))
display(modesPlot)

# Phugoid mode animation at highest payload
phugoidMode = 2
phugoidAnim = plot_mode_shapes_animation(eigenProblem[end],element2centralize=15,scale=20,nFramesPerCycle=41,fps=20,view=(30,30),showLegend=false,modeShapes=modeShapesOsc[phugoidMode:phugoidMode,end],plotSteady=false,plotBCs=false,plotAxes=false,displayProgress=true,save=true,savePath=string(relPath,"/heliosFlutterPRange_phugoid.gif"))
display(phugoidAnim)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    plot!(PRange, freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Payload [lb]", ylabel="Damping ratio", ylims=[-0.25,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(PRange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosFlutterPRange_Vgf.pdf"))

# Root locus 
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(dampsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosFlutterPRange_rootlocus.pdf"))

# Root locus (zoom)
plt3 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10])
for mode in 1:nModes
    scatter!(dampsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt3)
savefig(string(absPath,"/heliosFlutterPRange_rootlocuszoom.pdf"))

# Root locus (phugoid zoom)
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.1,0.2], ylims=[0,0.6])
for mode in 1:nModes
    scatter!(dampsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/heliosFlutterPRange_Phzoom.pdf"))

println("Finished heliosFlutterPRange.jl")