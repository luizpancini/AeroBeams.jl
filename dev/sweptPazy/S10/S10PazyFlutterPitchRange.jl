using AeroBeams, Plots, ColorSchemes

# Tip mass configurations (assumed tipMassPosOffset ahead/behind the LE/TE)
tipMassConfigs = ["LE", "TE"]
tipMasses = [11e-3, 11e-3]
tipMassPosOffsets = [0.05, 0.05]

# Root pitch angle range
θRange = π/180*[0,1,3,5,7,10]

# Airspeed range
URange = collect(1:1:100)

# Sweep angle [rad]
Λ = 10*π/180

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Number of modes
nModes = 5

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(tipMassConfigs),length(θRange),nModes)
tipOOP = Array{Float64}(undef,length(tipMassConfigs),length(θRange),length(URange))

# Sweep tip weight configurations
for (c,config) in enumerate(tipMassConfigs)
    # Set current tip mass position along the chord
    tipMassPos = config == "LE" ? (chord*normSparPos + tipMassPosOffsets[c]) : (- (chord*(1-normSparPos) + tipMassPosOffsets[c]))
    # Sweep pitch angle
    for (i,θ) in enumerate(θRange)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for for tip mass at $config, θ = $(round(θ*180/π)) deg, U = $U m/s")
            # Model
            model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMasses[c],ηtipMass=[0;tipMassPos;0])
            # Create and solve problem
            problem = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[2π,Inf])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[c,i,j] = problem.frequenciesOscillatory
            untrackedDamps[c,i,j] = problem.dampingsOscillatory
            untrackedEigenvectors[c,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP deflection
            tipOOP[c,i,j] = -problem.nodalStatesOverσ[end][nElem].u_n2[1]
        end
        # Apply mode tracking
        freqs[c,i,:],damps[c,i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[c,i,:],untrackedDamps[c,i,:],untrackedEigenvectors[c,i,:])
        # Separate frequencies and damping ratios by mode
        for mode in 1:nModes
            modeFrequencies[c,i,mode] = [freqs[c,i,j][mode] for j in eachindex(URange)]
            modeDampings[c,i,mode] = [damps[c,i,j][mode] for j in eachindex(URange)]
            modeDampingRatios[c,i,mode] = modeDampings[c,i,mode]./modeFrequencies[c,i,mode]
        end
    end
end

# Set paths
relPath = "/dev/sweptPazy/S10/outputs/S10PazyFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# V-g-f
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,105], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
for (c,config) in enumerate(tipMassConfigs)
    for (i,θ) in enumerate(θRange)
        for mode in 1:nModes
            plot!(URange, modeFrequencies[c,i,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        end
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for (i,θ) in enumerate(θRange)
        for mode in 1:nModes
            plot!(URange, modeDampingRatios[c,i,mode], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        end
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vf)
    display(plt_Vg)
    display(plt_Vgf)
    savefig(plt_Vf,string(absPath,"/S10PazyFlutterPitchRange_freq_",config,".pdf"))
    savefig(plt_Vg,string(absPath,"/S10PazyFlutterPitchRange_damp_",config,".pdf"))
    savefig(plt_Vgf,string(absPath,"/S10PazyFlutterPitchRange_Vgf_",config,".pdf"))
end

println("Finished S10PazyFlutterPitchRange.jl")