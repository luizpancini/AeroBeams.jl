using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,0.5:0.5:50)

# Flare angle [rad]
Λ = deg2rad(30)

# Pitch angle [rad]
θ = deg2rad(3)

# Sideslip angle [rad]
β = deg2rad(0)

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
eigensols = Array{EigenSolution}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U=$U m/s")
    # Update model
    model = create_HealySideslipFFWT(flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,flightDirection=[sin(β),cos(β),0])
    # Create and solve problem
    problem = create_EigenProblem(model=model,systemSolver=systemSolver,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e-2), maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTflutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 12
fs = 14
lfs = 12
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,165], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
for mode in 1:nModes
    plot!(URange, rads2Hz.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.1,0.05], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)

println("Finished HealySideslipFFWTflutterURange.jl")