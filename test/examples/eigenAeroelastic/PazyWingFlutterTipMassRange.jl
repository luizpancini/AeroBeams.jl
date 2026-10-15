using AeroBeams

# Get wing properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass configurations and respective positions
configurations = [1; 2; 3]
tipMassRange = [0; 0.01; 0.01]
tipMassPosRange = [0; chord*normSparPos; -chord*(1-normSparPos)]

# Airspeed range [m/s]
URange = vcat(1e-2,0.5:0.5:100)

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = true

# Pitch angle [rad]
θ = deg2rad(0)

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(configurations),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(configurations))
dampsOsc = Array{Matrix{Float64}}(undef,length(configurations))

# Sweep configurations
for (i,c) in enumerate(configurations)
    # Initial guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for configuration $c, U = $U m/s")
        # Model
        PazyWingFlutterTipMassRange,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,pitchAngle=θ,airspeed=U,tipMass=tipMassRange[i],ηtipMass=[0,tipMassPosRange[i],0])
        # Create and solve problem
        problem = create_EigenProblem(model=PazyWingFlutterTipMassRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(rads2Hz.(freqs),damps; fmin=1, fmax=60, maxTracks=nModes)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingFlutterTipMassRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
lstyles = [:solid :dash :dot]
gr()

# V-g-f for all configurations
plt11 = plot(ylabel="Frequency [Hz]")
for c in configurations
    for mode in 1:nModes
        plot!(URange, freqsOsc[c][mode,:], c=modeColors[mode], ls=lstyles[c], lw=lw, label=false)
    end
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.1,0.05], legend=:topleft)
plot!([NaN], [NaN], c=:black, ls=lstyles[1], lw=lw, label="Clean wing")
plot!([NaN], [NaN], c=:black, ls=lstyles[2], lw=lw, label="LE weight")
plot!([NaN], [NaN], c=:black, ls=lstyles[3], lw=lw, label="TE weight")
for c in configurations
    for mode in 1:nModes
        plot!(URange, dampsOsc[c][mode,:]./Hz2rads.(freqsOsc[c][mode,:]), c=modeColors[mode], ls=lstyles[c], lw=lw, label=false)
    end
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingFlutterTipMassRange.pdf"))

println("Finished PazyWingFlutterTipMassRange.jl")