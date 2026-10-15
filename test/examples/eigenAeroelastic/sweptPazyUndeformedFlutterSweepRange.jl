using AeroBeams

# Sweep angle range [rad]
ΛRange = deg2rad.([0,10,20,30])

# Airspeed range [m/s]
URange = collect(0:1:120)

# Root pitch angle [rad]
θ = deg2rad(0)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(ΛRange),length(URange))
tipOOP = Array{Float64}(undef,length(ΛRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(ΛRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(ΛRange))

# Sweep angle of sweep
for (i,Λ) in enumerate(ΛRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for Λ = $(round(Int,rad2deg(Λ))) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepAngle=Λ,pitchAngle=θ,airspeed=U,gravityVector=[0,0,-g],sweepStructuralCorrections=sweepStructuralCorrections)
        # Create and solve problem
        problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
        # Tip OOP deflection
        tipOOP[i,j] = problem.loadstepHistory.states[end][nElem].nodal.u2_b[3]
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(2), fmax=Hz2rads(50), maxTracks=nModes)
end

# Load reference data (from AePW4 meetings)
undef_damp_Lambda0_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_damp_Lambda0_NastranTechnion.txt")
undef_damp_Lambda20_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_damp_Lambda20_NastranTechnion.txt")
undef_damp_Lambda30_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_damp_Lambda30_NastranTechnion.txt")
undef_freq_Lambda0_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_freq_Lambda0_NastranTechnion.txt")
undef_freq_Lambda20_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_freq_Lambda20_NastranTechnion.txt")
undef_freq_Lambda30_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/undef_freq_Lambda30_NastranTechnion.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazyUndeformedFlutterSweepRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad([:blue,:red,:green], nModes, categorical=true)
ts = 10
fs = 15
lfs = 8
lw = 2
ms = 4
gr()

# V-g-f
for (i,Λ) in enumerate(ΛRange)
    Λstr = string(round(Int,rad2deg(Λ)))
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,120], ylims=[0,50], xticks=0:10:120, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legendposition=(0.75,0.95))
    if i==1
        plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
        scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Nastran (Technion)")
    end
    for mode in 1:nModes
        plot!(URange, rads2Hz(freqsOsc[i][mode,:]), c=modeColors[mode], lw=lw, label=false)
    end
    if i == 1
        scatter!(undef_freq_Lambda0_NastranTechnion[1,1:46],undef_freq_Lambda0_NastranTechnion[2,1:46], c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda0_NastranTechnion[1,48:93],undef_freq_Lambda0_NastranTechnion[2,48:93], c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda0_NastranTechnion[1,95:end],undef_freq_Lambda0_NastranTechnion[2,95:end], c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 3
        scatter!(undef_freq_Lambda20_NastranTechnion[1,1:31],undef_freq_Lambda20_NastranTechnion[2,1:31], c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda20_NastranTechnion[1,33:63],undef_freq_Lambda20_NastranTechnion[2,33:63], c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_freq_Lambda20_NastranTechnion[1,65:end],undef_freq_Lambda20_NastranTechnion[2,65:end], c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 4
        scatter!(undef_freq_Lambda30_NastranTechnion[1,1:31],undef_freq_Lambda30_NastranTechnion[2,1:31], c=modeColors[1], ms=ms, msw=msw, label=false)    
        scatter!(undef_freq_Lambda30_NastranTechnion[1,33:63],undef_freq_Lambda30_NastranTechnion[2,33:63], c=modeColors[2], ms=ms, msw=msw, label=false)  
        scatter!(undef_freq_Lambda30_NastranTechnion[1,65:end],undef_freq_Lambda30_NastranTechnion[2,65:end], c=modeColors[3], ms=ms, msw=msw, label=false)  
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,120], ylims=[-0.15,0.075], xticks=0:10:120, yticks=-0.1:0.05:0.1, tickfont=font(ts), guidefont=font(fs))
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URange, dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=modeColors[mode], lw=lw, label=false)
    end
    if i == 1
        scatter!(undef_damp_Lambda0_NastranTechnion[1,1:6],undef_damp_Lambda0_NastranTechnion[2,1:6]/2, c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda0_NastranTechnion[1,8:37],undef_damp_Lambda0_NastranTechnion[2,8:37]/2, c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda0_NastranTechnion[1,39:end],undef_damp_Lambda0_NastranTechnion[2,39:end]/2, c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 3
        scatter!(undef_damp_Lambda20_NastranTechnion[1,1:4],undef_damp_Lambda20_NastranTechnion[2,1:4]/2, c=modeColors[1], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda20_NastranTechnion[1,6:33],undef_damp_Lambda20_NastranTechnion[2,6:33]/2, c=modeColors[2], ms=ms, msw=msw, label=false)
        scatter!(undef_damp_Lambda20_NastranTechnion[1,35:end],undef_damp_Lambda20_NastranTechnion[2,35:end]/2, c=modeColors[3], ms=ms, msw=msw, label=false)
    elseif i == 4
        scatter!(undef_damp_Lambda30_NastranTechnion[1,1:5],undef_damp_Lambda30_NastranTechnion[2,1:5]/2, c=modeColors[1], ms=ms, msw=msw, label=false)    
        scatter!(undef_damp_Lambda30_NastranTechnion[1,6:36],undef_damp_Lambda30_NastranTechnion[2,6:36]/2, c=modeColors[2], ms=ms, msw=msw, label=false) 
        scatter!(undef_damp_Lambda30_NastranTechnion[1,38:end],undef_damp_Lambda30_NastranTechnion[2,38:end]/2, c=modeColors[3], ms=ms, msw=msw, label=false) 
    end
    annotate!([30],[0.05], text("\$\\Lambda=$Λstr ^\\circ\$", 30, :black))
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(string(absPath,"/sweptPazyUndeformedFlutterSweepRange_Vgf_Lambda",Λstr,".pdf"))
end

println("Finished sweptPazyUndeformedFlutterSweepRange.jl")