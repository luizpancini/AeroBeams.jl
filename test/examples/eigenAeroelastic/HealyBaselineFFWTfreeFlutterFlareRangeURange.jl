using AeroBeams

# Hinge configuration
hingeConfiguration = :free

# Flare angle range [rad]
ΛRange = deg2rad.([10,15,20])

# Airspeed range [m/s]
URange = vcat(0.1,0.5:0.5:40)

# Pitch angle [rad]
θ = deg2rad(0)

# Gravity
g = 0

# Discretization
nElementsInner = 16
nElementsFFWT = 4

# Tip loss options (assumed, since Healy's analysis uses DLM for aerodynamic)
hasTipCorrection = true
tipLossDecayFactor = 12

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Number of modes to track
nModes = 2

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(ΛRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(ΛRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(ΛRange))

# Sweep flare angle
for (i,Λ) in enumerate(ΛRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for Λ=$(round(Int,rad2deg(Λ))) deg, U=$U m/s")
        # Update model
        model = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,g=g,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
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
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e-2), maxTracks=nModes)
end

# Load reference data
flare10_mode1_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode1_damp.txt")
flare10_mode2_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode2_damp.txt")
flare10_mode1_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode1_freq.txt")
flare10_mode2_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode2_freq.txt")

flare15_mode1_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode1_damp.txt")
flare15_mode2_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode2_damp.txt")
flare15_mode1_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode1_freq.txt")
flare15_mode2_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode2_freq.txt")

flare20_mode1_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode1_damp.txt")
flare20_mode2_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode2_damp.txt")
flare20_mode1_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode1_freq.txt")
flare20_mode2_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode2_freq.txt")

healy = [
    (f1=flare10_mode1_freq, f2=flare10_mode2_freq, d1=flare10_mode1_damp, d2=flare10_mode2_damp),
    (f1=flare15_mode1_freq, f2=flare15_mode2_freq, d1=flare15_mode1_damp, d2=flare15_mode2_damp),
    (f1=flare20_mode1_freq, f2=flare20_mode2_freq, d1=flare20_mode1_damp, d2=flare20_mode2_damp),
]

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTfreeFlutterFlareRangeURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
lw = 2

# V-g-f
plt1 = plot(ylabel="Frequency [Hz]", xlims=[0,40], ylims=[0,5], legend=:bottomright)
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,40], ylims=[-0.5,0.25], legend=:bottomleft)
plot!(plt1, [NaN],[NaN], lc=:black, ls=:dash,  lw=lw, label="Healy (2023)")
plot!(plt1, [NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,Λ) in enumerate(ΛRange)
    for mode in 1:nModes
        plot!(plt1, URange, rads2Hz.(freqsOsc[i][mode,:]), c=colors[i], lw=lw, label=false)
        plot!(plt2, URange, dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=colors[i], lw=lw, label=false)
    end
    h = healy[i]
    plot!(plt1, h.f1[1,:], h.f1[2,:], c=colors[i], ls=:dash, lw=lw, label=false)
    plot!(plt1, h.f2[1,:], h.f2[2,:], c=colors[i], ls=:dash, lw=lw, label=false)
    plot!(plt2, h.d1[1,:], h.d1[2,:], c=colors[i], ls=:dash, lw=lw, label=false)
    plot!(plt2, h.d2[1,:], h.d2[2,:], c=colors[i], ls=:dash, lw=lw, label=false)
    plot!(plt2, [NaN],[NaN], c=colors[i], ls=:solid, lw=lw, label="\$\\Lambda = $(round(Int,rad2deg(ΛRange[i]))) \\degree\$")
end
plt = plot(plt1, plt2, layout=(2,1))
display(plt)

println("Finished HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")