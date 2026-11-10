using AeroBeams

# Tip mass configurations
tipMassConfigs = [:LE, :TE]

# Root pitch angle range [rad]
θRange = deg2rad.([0,1,3,5,7,8.5,10])

# Airspeed range [m/s]
URange = collect(30:1:85)

# Sweep angle [rad]
Λ = deg2rad(20)

# Flag to model the wingtip rod (true) or NASTRAN with stiffened wingtip rod (false)
flexibleWingtipRod = true

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = flexibleWingtipRod

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(Λ=Λ)

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(tipMassConfigs),length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(tipMassConfigs),length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs),length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs),length(θRange))
flutterOnsetSpeed = fill(Inf,length(tipMassConfigs),length(θRange))
flutterOnsetFreq = fill(NaN,length(tipMassConfigs),length(θRange))
flutterOnsetTipOOP = fill(NaN,length(tipMassConfigs),length(θRange))

# Sweep tip mass configurations
for (c,config) in enumerate(tipMassConfigs)
    # Sweep root pitch angle
    for (i,θ) in enumerate(θRange)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for tip mass at $config, θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
            # Model
            model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepAngle=Λ,pitchAngle=θ,airspeed=U,sweptTipMassConfig=config,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)
            # Create and solve problem
            problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,x0=x0)
            solve!(problem)
            # Eigensolution
            eigensols[c,i,j] = problem.eigensolution
            # Tip OOP deflection
            tipOOP[c,i,j] = -problem.loadstepHistory.states[end][nElem].nodal.u2[1]
            # Set initial guess solution for next iteration
            x0 = problem.systemArrays.x
        end
        # Apply mode tracking
        freqs,damps,_ = mode_tracking(eigensols[c,i,:])
        # Get frequencies and dampings of oscillatory modes
        freqsOsc[c,i],dampsOsc[c,i],_ = filter_tracked_modes(freqs,damps; fmin=2π, maxTracks=nModes)
        # Flutter onset speed, frequency and tip OOP (lowest onset)
        onsetSpeeds,onsetFreqs,_,onsetTipOOPs,_ = find_flutter(URange,freqsOsc[c,i],dampsOsc[c,i],tipOOP[c,i,:])
        if !isempty(onsetSpeeds)
            flutterOnsetSpeed[c,i] = onsetSpeeds[1]
            flutterOnsetFreq[c,i] = onsetFreqs[1]
            flutterOnsetTipOOP[c,i] = onsetTipOOPs[1]
        end
    end
end

# Load reference data (experiments from Technion - https://doi.org/10.5281/zenodo.16354530)
Uf_vs_tipdisp_LE_c = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_LE_c.txt")
Uf_vs_tipdisp_LE_g = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_LE_g.txt")
Uf_vs_tipdisp_TE_c = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_TE_c.txt")
Uf_vs_tipdisp_TE_g = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Uf_vs_tipdisp_TE_g.txt")
Ff_vs_tipdisp_LE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Ff_vs_tipdisp_LE.txt")
Ff_vs_tipdisp_TE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/Ff_vs_tipdisp_TE.txt")

## Plots
# Set paths
relPath = "/dev/Pazy/S20/outputs/S20PazyFlutterOnsetPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
configColors = cgrad(:rainbow, length(tipMassConfigs), categorical=true)
ts = 10
fs = 10
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# Flutter onset speed vs. tip displacement
plt1 = plot(ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,80], xticks=vcat(0:5:40), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:black, label="Exp.")
plot!([NaN], [NaN], lw=lw, c=configColors[1], label="LE weight")
plot!([NaN], [NaN], lw=lw, c=configColors[2], label="TE weight")
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:]/L*100, flutterOnsetSpeed[c,:], c=configColors[c], lw=lw, label=false)
    if config == :LE
        plot!(Uf_vs_tipdisp_LE_g[1,:]/L*100, Uf_vs_tipdisp_LE_g[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == :TE
        plot!(Uf_vs_tipdisp_TE_g[1,:]/L*100, Uf_vs_tipdisp_TE_g[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end

# Flutter onset frequency vs. tip displacement
plt2 = plot(xlabel="Normalized tip displacement [% semispan]", ylabel="Flutter frequency [Hz]", xlims=[0,40], ylims=[0,35], xticks=vcat(0:5:40), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:]/L*100, rads2Hz.(flutterOnsetFreq[c,:]), c=configColors[c], lw=lw, label=false)
    if config == :LE
        plot!(Ff_vs_tipdisp_LE[1,:]*100, Ff_vs_tipdisp_LE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == :TE
        plot!(Ff_vs_tipdisp_TE[1,:]*100, Ff_vs_tipdisp_TE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end

# Combine plots and save
plt = plot(plt1, plt2, layout=(2,1))
display(plt)
savefig(plt,string(absPath,"/S20PazyFlutterOnsetPitchRange_Uf_omegaf.pdf"))

println("Finished S20PazyFlutterOnsetPitchRange.jl")