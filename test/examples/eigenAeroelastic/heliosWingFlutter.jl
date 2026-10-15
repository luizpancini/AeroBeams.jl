using AeroBeams

# Airspeed range [m/s]
URange = collect(0.1:0.1:15)

# Stiffness factor
λ = 1

# Flag to include beam pod
beamPods = false

# Aerodynamic solver
aeroSolver = Indicial()

# Number of modes to track
nModes = 8

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
problem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))

# Initial guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model 
    _,_,heliosWing,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,stiffnessFactor=λ,airspeed=U)
    # Set initial guess solution
    x0 = i==1 ? zeros(0) : problem[i-1].systemArrays.x
    # Create and solve trim problem
    problem[i] = create_EigenProblem(model=heliosWing,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, fmax=20, maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/heliosWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at the last airspeed
modesPlot = plot_mode_shapes(problem[end],modeShapes=modeShapesOsc[:,end],modeNumbers=collect(1:nModes),scale=5,view=(30,30),save=true,savePath=string(relPath,"/heliosWingFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
msw = 0
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]", ylims=[0,15])
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-2.0,0.25], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosWingFlutter_Vgf.pdf"))

# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-0.5,0.25] ,ylims=[0,15])
for mode in 1:nModes
    plot!(dampsOsc[mode,:]./freqsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosWingFlutter_rootlocus.pdf"))

println("Finished heliosWingFlutter.jl")