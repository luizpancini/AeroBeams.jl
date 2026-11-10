using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:130)

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = true

# Root pitch angle [rad]
θ = deg2rad(7)

# Mesh: :nominal or :refined
mesh = :refined

# Initialize model
PazyWingNoRodFlutter,nElem,L,chord,normSparPos = create_Pazy(aeroSolver=aeroSolver,withWingtipRod=false,airfoil=airfoil,upright=upright,pitchAngle=θ,mesh=mesh)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = if aeroSolver isa DynamicStallSolver
    mesh == :nominal ? 150 : 200
else
    mesh == :nominal ? 45 : 65
end

# Initialize outputs
problem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update airspeed
    update_v_A!(PazyWingNoRodFlutter,v_A=[0,U,0])
    # Set initial guess solution
    x0 = i==1 ? zeros(0) : problem[i-1].systemArrays.x
    # Create and solve problem
    problem[i] = create_EigenProblem(model=PazyWingNoRodFlutter,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem[i])    
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)
rads2Hz!(freqs)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e0, maxTracks=nModes)

# Flutter onset and offset speeds, frequencies, modes and tip OOP displacements
flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetModes,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetModes = find_flutter(URange,freqsOsc,dampsOsc)

# Show flutter results
for (m,speed,freq) in zip(flutterOnsetModes,flutterOnsetSpeeds,flutterOnsetFreqs)
    println("Mode $m: onset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) Hz")
end
for (m,speed,freq) in zip(flutterOffsetModes,flutterOffsetSpeeds,flutterOffsetFreqs)
    println("Mode $m: offset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) Hz")
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingNoRodFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]", ylims=[0,130])
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.1], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./Hz2rads.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingNoRodFlutter_Vgf.pdf"))

println("Finished PazyWingNoRodFlutter.jl")