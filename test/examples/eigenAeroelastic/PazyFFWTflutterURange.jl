using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:40)

# Hinge node
hingeNode = 13

# Flare angle [rad]
Λ = deg2rad(10)

# Pitch angle [rad]
θ = deg2rad(3)

# Sideslip angle [rad]
β = deg2rad(0)

# Tip mass [kg] and its position [m]
tipMass = 10e-3
ηtipMass = 1e-3*[0,50,0]

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
ϕHinge = Array{Float64}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U=$U m/s")
    # Update model
    model,_ = create_Pazy(hinged=true,hingeNode=hingeNode,flareAngle=Λ,airspeed=U,pitchAngle=θ,flightDirection=[-sin(β),cos(β),0],tipMass=tipMass,ηtipMass=ηtipMass)
    # Create and solve problem
    problem = create_EigenProblem(model=model,systemSolver=systemSolver,x0=x0,eigenmethod=eigenmethod,nev=nev)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Coast angle
    ϕHinge[i] = rad2deg(problem.loadstepHistory.hacStates[end][1].ϕ)
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e-2), maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyFFWTflutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2

# V-g-f
plt1 = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,120], legend=:bottomright)
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.1,0.05], legend=:bottomleft)
for m in 1:nModes
    plot!(plt1, URange, rads2Hz.(freqsOsc[m,:]), c=colors[m], lw=lw, label=false)
    plot!(plt2, URange, dampsOsc[m,:]./freqsOsc[m,:], c=colors[m], lw=lw, label=false)
end
plt = plot(plt1, plt2, layout=(2,1))
display(plt)
savefig(string(absPath,"/PazyFFWTflutterURange_Vgf.pdf"))

println("Finished PazyFFWTflutterURange.jl")