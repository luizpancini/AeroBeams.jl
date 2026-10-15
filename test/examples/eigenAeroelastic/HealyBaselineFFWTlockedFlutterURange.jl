using AeroBeams

# Hinge configuration
hingeConfiguration = :locked

# Airspeed range [m/s]
URange = collect(0:1:120)

# Pitch angle [rad]
θ = deg2rad(7)

# Discretization
nElementsInner = 16
nElementsFFWT = 4

# Tip loss options (assumed, since Healy's analysis uses DLM for aerodynamic)
hasTipCorrection = true
tipLossDecayFactor = 12

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U=$U m/s")
    # Update model
    model =  create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
    # Create and solve problem
    problem = create_EigenProblem(model=model,systemSolver=systemSolver,eigenmethod=eigenmethod,nev=nev)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e-1), maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTlockedFlutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2

# V-g-f
plt1 = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,50], legend=:bottomright)
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.2,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(plt1, URange, rads2Hz.(freqsOsc[mode,:]), c=colors[mode], lw=lw, label=false)
    plot!(plt2, URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=colors[mode], lw=lw, label=false)
end
plt = plot(plt1, plt2, layout=(2,1))
display(plt)

println("Finished HealyBaselineFFWTlockedFlutterURange.jl")