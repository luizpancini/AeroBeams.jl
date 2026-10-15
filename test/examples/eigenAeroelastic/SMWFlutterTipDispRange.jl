using AeroBeams

# Tip force range [N]
F3Range = collect(-50:2:50)

# Airspeed range [m/s]
URange = collect(20:0.5:35)

# Altitude [m]
h = 20e3

# Pitch angle [rad]
θ = deg2rad(0)

# Gravity
g = 0

# Discretization
nElem = 16

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil
airfoil = deepcopy(FlatPlate)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(F3Range),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(F3Range))
dampsOsc = Array{Matrix{Float64}}(undef,length(F3Range))
flutterSpeed = Array{Float64}(undef,length(F3Range))
flutterFreq = Array{Float64}(undef,length(F3Range))
flutterTipDisp = Array{Float64}(undef,length(F3Range))

# Initialize guess solution
x0 = zeros(0)

# Sweep tip force
for (j,F3) in enumerate(F3Range)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (k,U) in enumerate(URange)
        # Display progress
        println("Solving F3 = $F3 N, U = $U m/s")
        # Update model
        SMWFlutterTipDispRange,_ = create_SMW(aeroSolver=aeroSolver,airfoil=airfoil,pitchAngle=θ,altitude=h,g=g,tipF3=F3,airspeed=U,nElem=nElem)
        # Create and solve problem
        problem = create_EigenProblem(model=SMWFlutterTipDispRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[k] = problem.eigensolution
        # Tip OOP displacement
        tip_u3[j,k] = problem.loadstepHistory.states[end][nElem].nodal.u2_b[3]
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols)
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[j],dampsOsc[j],_ = filter_tracked_modes(freqs,damps; fmin=1e-2, maxTracks=nModes, tolDamp=1e-9)
    # Flutter speeds, frequencies and tip displacements
    flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterOnsetTipDisp,_ = find_flutter(URange,freqsOsc[j],dampsOsc[j],tip_u3[j,:])
    # Set flutter speed and corresponding flutter frequency and tip displacement
    flutterSpeed[j] = flutterOnsetSpeeds[1]
    flutterFreq[j] = flutterOnsetFreqs[1]
    flutterTipDisp[j] = flutterOnsetTipDisp[1]
end

# Load reference data
flutterSpeedRef = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipDisp.txt"))
flutterFreqRef = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterFreqVsTipDisp.txt"))
flutterSpeedVsDispRef = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipDispFull.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWFlutterTipDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 3
gr()

# Flutter speed and frequency vs. tip displacement
plt11 = plot(ylabel="Flutter speed [m/s]", xlims=[0,3], ylims=[0,35], legend=:bottomleft)
plot!(flutterTipDisp, flutterSpeed, c=:black, lw=lw, label="AeroBeams")
plot!(flutterSpeedRef[1,:], flutterSpeedRef[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plt12 = plot(xlabel="Tip displacement [m]", ylabel="Flutter frequency [rad/s]", xlims=[0,3], ylims=[0,35])
plot!(flutterTipDisp, flutterFreq, c=:black, lw=lw, label=false)
plot!(flutterFreqRef[1,:], flutterFreqRef[2,:], c=:black, ls=:dash, lw=lw, label=false)
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWFlutterTipDispRange_flutterVsDisp.pdf"))

# Complete flutter speed curve
plt2 = plot(xlabel="Tip displacement [m]", ylabel="Flutter speed [m/s]", xlims=[-3,3], ylims=[0,35], legend=:bottomleft)
plot!(flutterTipDisp, flutterSpeed, c=:black, lw=lw, label="AeroBeams")
plot!(flutterSpeedVsDispRef[1,:], flutterSpeedVsDispRef[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
display(plt2)
savefig(string(absPath,"/SMWFlutterTipDispRange_flutterCurve.pdf"))

println("Finished SMWFlutterTipDispRange.jl")