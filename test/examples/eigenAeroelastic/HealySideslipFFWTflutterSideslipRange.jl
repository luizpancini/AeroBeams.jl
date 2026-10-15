using AeroBeams

# Pitch angle range [rad]
θRange = deg2rad.([-3,3,9,15])

# Sideslip angle range [rad]
βRange = deg2rad.(vcat(-20:1:20))

# Flare angle [rad]
Λ = deg2rad(30)

# Airspeed [m/s]
U = 25

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 10

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(θRange),length(βRange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))

# Sweep root pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep sideslip angle
    for (j,β) in enumerate(βRange)
        # Display progress
        println("Solving for θ=$(round(Int,rad2deg(θ))) deg, β=$(round(Int,rad2deg(β))) deg")
        # Update model
        model = create_HealySideslipFFWT(flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,flightDirection=[-sin(β),cos(β),0])
        # Create and solve problem
        problem = create_EigenProblem(model=model,systemSolver=systemSolver,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
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
relPath = "/test/outputs/figures/HealySideslipFFWTflutterSideslipRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2

# V-g-f
plt1 = plot(ylabel="Frequency [Hz]", xlims=extrema(rad2deg.(βRange)), ylims=[0,50], legend=:bottomright)
plt2 = plot(xlabel="Sideslip Angle [deg]", ylabel="Damping Ratio", xlims=extrema(rad2deg.(βRange)), ylims=[-0.1,0.05], legend=:bottomleft)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(plt1, rad2deg.(βRange), rads2Hz.(freqsOsc[i][mode,:]), c=colors[i], lw=lw, label=false)
        plot!(plt2, rad2deg.(βRange), dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=colors[i], lw=lw, label=false)
    end
end
plt = plot(plt1, plt2, layout=(2,1))
display(plt)

println("Finished HealySideslipFFWTflutterSideslipRange.jl")