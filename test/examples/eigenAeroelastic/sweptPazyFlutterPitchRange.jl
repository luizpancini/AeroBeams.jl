using AeroBeams

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass configuration (select between :none, :TE or :LE)
tipMassConfig = :TE

# Corresponding tip mass and its position over the chord
if tipMassConfig == :none
    tipMass = 0
    tipMassPos = 0
elseif tipMassConfig == :TE
    tipMass = 15.5e-3
    tipMassPos = -chord*(1-normSparPos) - 0.01 # Assume about 1 cm behind TE
elseif tipMassConfig == :LE
    tipMass = 15.5e-3
    tipMassPos = chord*normSparPos + 0.01 # Assume about 1 cm ahead of LE
else
    error("Select tipMassConfig as :none, :TE or :LE")
end

# Sweep angle [rad]
Λ = deg2rad(20)

# Root pitch angle range [rad]
θRange = deg2rad.([0,1,3,5,7,10])

# Airspeed range [m/s]
URange = collect(0:0.5:100)

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
upright = true

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(Int,rad2deg(θ))) deg, U = $U m/s")
        # Model
        sweptPazyFlutterPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepAngle=Λ,pitchAngle=θ,airspeed=U,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0])
        # Create and solve problem
        problem = create_EigenProblem(model=sweptPazyFlutterPitchRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
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
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=2π, maxTracks=nModes)
end

# Load reference data (experiments from Technion)
freqs_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa0.txt")
freqs_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa1.txt")
freqs_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa3.txt")
freqs_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa5.txt")
freqs_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa7.txt")
freqs_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa10.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazyFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
tipLossStr = ifelse(hasTipCorrection, tipLossType, "none")
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# V-g-f
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,45], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
scatter!([NaN], [NaN], ms=ms, msw=msw, mc=:white, msc=:black, label="Exp.")
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, rads2Hz.(freqsOsc[i][mode,:]), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        if i==1 && mode <= size(freqs_aoa0_ref,1)-1
            scatter!(freqs_aoa0_ref[1,:], freqs_aoa0_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif i==2 && mode <= size(freqs_aoa1_ref,1)-1
            scatter!(freqs_aoa1_ref[1,:], freqs_aoa1_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif i==3 && mode <= size(freqs_aoa3_ref,1)-1
            scatter!(freqs_aoa3_ref[1,:], freqs_aoa3_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif i==4 && mode <= size(freqs_aoa5_ref,1)-1
            scatter!(freqs_aoa5_ref[1,:], freqs_aoa5_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif i==5 && mode <= size(freqs_aoa7_ref,1)-1
            scatter!(freqs_aoa7_ref[1,:], freqs_aoa7_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif i==6 && mode <= size(freqs_aoa10_ref,1)-1
            scatter!(freqs_aoa10_ref[1,:], freqs_aoa10_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)      
        end
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vf)
display(plt_Vg)
display(plt_Vgf)
savefig(plt_Vf,string(absPath,"/sweptPazyFlutterPitchRange_freq_",tipLossStr,"_Lambda",string(round(Int,rad2deg(Λ))),".pdf"))
savefig(plt_Vg,string(absPath,"/sweptPazyFlutterPitchRange_damp_",tipLossStr,"_Lambda",string(round(Int,rad2deg(Λ))),".pdf"))
savefig(plt_Vgf,string(absPath,"/sweptPazyFlutterPitchRange_Vgf_",tipLossStr,"_Lambda",string(round(Int,rad2deg(Λ))),".pdf"))

println("Finished sweptPazyFlutterPitchRange.jl")