using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,0.1:0.1:40)

# Altitude
h = 20e3

# Pitch angle
θ = 0

# Discretization
nElem = 16

# Aerodynamic solver
aeroSolver = Indicial()

# Model
SMWLinearFlutter,L = create_SMW(aeroSolver=aeroSolver,altitude=h,pitchAngle=θ,nElem=nElem)

# Linear system solver
systemSolver = LinearSolver()

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 30

# Initialize outputs
problem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update velocity of basis A 
    update_v_A!(SMWLinearFlutter,v_A=[0,U,0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=SMWLinearFlutter,systemSolver=systemSolver,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory and non-oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=eps(), fmax=1e2, maxTracks=nModes)
freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps; fmin =0, fmax=0)

# Flutter speed, frequency and modes
flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterInds = find_flutter(URange,freqsOsc,dampsOsc,1:length(URange))
flutterSpeed = minimum(flutterOnsetSpeeds)
flutterFreq = flutterOnsetFreqs[argmin(flutterOnsetSpeeds)]

## Divergence speed
divergenceSpeeds,_ = find_divergence(URange,freqsNonOsc,dampsNonOsc)
divergenceSpeed = minimum(divergenceSpeeds)

# Reference solution by Patil & Hodges & Cesnik: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)
flutterSpeedRef = 32.21
flutterFreqRef = 22.61
divergenceSpeedRef = 37.29

# Compute relative differences [%]
ϵUf = 100*(flutterSpeed/flutterSpeedRef - 1)
ϵFf = 100*(flutterFreq/flutterFreqRef - 1)
ϵUd = 100*(divergenceSpeed/divergenceSpeedRef - 1)

# Pretty table of results
data = [
    "Divergence speed"   divergenceSpeed  divergenceSpeedRef ϵUd
    "Flutter speed"      flutterSpeed  flutterSpeedRef  ϵUf
    "Flutter frequency"  flutterFreq   flutterFreqRef   ϵFf
]
pretty_table(
    data;
    column_labels = ["", "AeroBeams", "Patil & Hodges & Cesnik", "Error [%]"],
    title         = "16-meter wing",
    formatters    = [
        (v, i, j) -> (v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j == 4 && v isa Number) ? string(round(v, sigdigits=2), " %") : v,
    ],
    source_notes  = "Reference: Patil & Hodges & Cesnik, Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)",
)

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWLinearFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes at flutter condition
iFlutter = round(Int,flutterInds[1])+1
modesPlot = plot_mode_shapes(problem[iFlutter],modeShapes=modeShapesOsc[1:nModes,iFlutter],modeNumbers=collect(1:nModes),scale=2,view=(30,30),legendPos=:best,save=true,savePath=string(relPath,"/SMWLinearFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 4
msw = 0
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    scatter!(URange, freqsOsc[mode,:], mc=modeColors[mode], ms=ms, msw=msw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.2,0.1])
for mode in 1:nModes
    scatter!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], mc=modeColors[mode], ms=ms, msw=msw, label=false)
    scatter!([NaN], [NaN], mc=modeColors[mode], ms=ms, msw=msw, label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWLinearFlutter_Vgf.pdf"))

println("Finished SMWLinearFlutter.jl")