using AeroBeams

# Wingtip configuration
wingtipTab = false

# Pitch angle range [rad]
θRange = deg2rad.([2.5, 7.5])

# Airspeed range [m/s]
URange = vcat(1e-3,0.5:0.5:30)

# Discretization
nElementsInner = 16
nElementsOuter = 6
nElementsFFWT = 8

# Tip loss flag
hasTipCorrection = true

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Number of modes of intereset
nModes = 8

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))
ϕHinge = Array{Float64}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress 
        println("Solving for θ=$(round(rad2deg(θ),digits=1)) deg, U=$U m/s") 
        # Update model
        model = create_HealyLCOFFWT(wingtipTab=wingtipTab,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)
        # Create and solve problem
        problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,systemSolver=systemSolver,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
        # Coast angle
        ϕHinge[i,j] = -rad2deg(problem.loadstepHistory.hacStates[end][1].ϕ)
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e-1), maxTracks=nModes)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyLCOFFWTfreeFlutterAoARange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2

# V-g-f
plt1 = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,60], legend=:bottomright)
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.25,0.05], legend=:bottomleft)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        label = mode==1 ? "\$\\theta=$(round(rad2deg(θ),digits=1))\$ deg" : false
        plot!(plt1, URange, rads2Hz.(freqsOsc[i][mode,:]), c=colors[i], lw=lw, label=false)
        plot!(plt2, URange, dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=colors[i], lw=lw, label=label)
    end
end
plt = plot(plt1, plt2, layout=(2,1))
display(plt)

println("Finished HealyLCOFFWTfreeFlutterAoARange.jl")