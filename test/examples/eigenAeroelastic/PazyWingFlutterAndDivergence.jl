using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,0.5:0.5:120)

# Pitch angle [rad]
θ = deg2rad(0)

# Gravity
g = 0

# Aerodynamic solver
aeroSolver = Inflow()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = true

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass
tipMass = 0.0
ηtipMass = [0; -chord*(1-normSparPos)-0.0; 0]

# Initialize model
PazyWingFlutterAndDivergence,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,pitchAngle=θ,gravityVector=[0,0,-g],tipMass=tipMass,ηtipMass=ηtipMass)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 90

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
tip_OOP = Array{Float64}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update velocity of basis A (and update model)
    update_v_A!(PazyWingFlutterAndDivergence,v_A=[0,U,0])
    # Create and solve problem
    problem = create_EigenProblem(model=PazyWingFlutterAndDivergence,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem)
    # Tip OOP displacement at midchord
    states = problem.loadstepHistory.states
    tipTwist = twist_angle(states[end][nElem].nodal.p2_b, outputUnit=:deg)
    tip_OOP[i] = -(states[end][nElem].nodal.u2[1] - chord*(1/2-normSparPos)*sind(tipTwist))/L*100
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)
rads2Hz!(freqs)

# Get frequencies and dampings of oscillatory and non-oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=1e-2, fmax=120, maxTracks=nModes, whichTracks=:last)
freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps; fmin=0, fmax=0)

# Divergence speed
divergenceSpeeds,_ = find_divergence(URange,freqsNonOsc,dampsNonOsc)
divergenceSpeed = minimum(divergenceSpeeds)
println("Divergence speed = $(round(divergenceSpeed, digits=2)) m/s")

# Flutter onset and offset speeds, frequencies, modes and tip OOP displacements
flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetModes,flutterOnsetDisp,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetModes,flutterOffsetDisp = find_flutter(URange,freqsOsc,dampsOsc,tip_OOP)

# Show flutter results
for (m,speed,freq,disp) in zip(flutterOnsetModes,flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp)
    println("Mode $m: onset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) Hz, tip OOP = $(round(disp, digits=2)) %")
end
for (m,speed,freq,disp) in zip(flutterOffsetModes,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp)
    println("Mode $m: offset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) Hz, tip OOP = $(round(disp, digits=2)) %")
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingFlutterAndDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]", xlims=[0,120], ylims=[0,50])
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,120], ylims=[-0.25,0.1], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./Hz2rads.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingFlutterAndDivergence_Vgf.pdf"))

println("Finished PazyWingFlutterAndDivergence.jl")