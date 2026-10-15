using AeroBeams, Plots, ColorSchemes

# Tip mass configuration and tip mass
tipMassConfig = "LE"
tipMass = 11e-3

# Tip mass position offset range
tipMassPosOffsetRange = 1e-2*collect(0:1:10)

# Airspeed range
URange = collect(1:1:100)

# Root pitch angle
θRange = 10π/180

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
nModes = 3

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(tipMassPosOffsetRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(tipMassPosOffsetRange),nModes)
tipOOP = Array{Float64}(undef,length(tipMassPosOffsetRange),length(URange))

# Sweep tip mass position offset range
for (i,o) in enumerate(tipMassPosOffsetRange)
    tipMassPos = tipMassConfig == "LE" ? (chord*normSparPos + tipMassPosOffset) : (- (chord*(1-normSparPos) + tipMassPosOffset))
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for for tip mass position offset $(o*1e2) cm, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;o;0])
        # Create and solve problem
        problem = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[2π,Inf])
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = problem.dampingsOscillatory
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        # Tip OOP deflection
        tipOOP[i,j] = -problem.nodalStatesOverσ[end][nElem].u_n2[1]
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

# Case ID
idString = string(tipMassConfig, "_mass", round(Int,tipMass*1e3), "_Theta", round(θ*180/π,digits=1))

# Set paths
relPath = "/dev/sweptPazy/S10/outputs/S10PazyFlutterTipMassPosRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
tipMassPosColors = cgrad(:rainbow, length(tipMassPosOffsetRange), categorical=true)
ts = 10
fs = 16
lw = 2
ms = 4
msw = 1
gr()

# V-g-f with selected modes
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,o) in enumerate(tipMassPosOffsetRange)
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(tipMassPosOffsetRange), label=false)
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,o) in enumerate(tipMassPosOffsetRange)
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[i,mode], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(tipMassPosOffsetRange), label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vf)
display(plt_Vg)
display(plt_Vgf)
savefig(plt_Vf,string(absPath,"/S10PazyFlutterTipMassPosRange_freq_",idString,".pdf"))
savefig(plt_Vg,string(absPath,"/S10PazyFlutterTipMassPosRange_damp_",idString,".pdf"))
savefig(plt_Vgf,string(absPath,"/S10PazyFlutterTipMassPosRange_Vgf_",idString,".pdf"))

# V-g of flutter mode with tip mass position offset as variable
flutterMode = 3
plt_Vg2 = plot(title="Tip mass @$tipMassConfig, mass $(tipMass*1e3) g, \$\\theta = $(round(θ*180/π,digits=1))\$ deg", xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[25,50], ylims=[-0.05,0.025], xticks=vcat(0:5:100), tickfont=font(ts), guidefont=font(fs), colorbar_title="Tip mass position [cm]")
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,m) in enumerate(tipMassPosOffsetRange*1e2)
    plot!(URange, modeDampingRatios[i,flutterMode], lz=m, c=tipMassPosColors, lw=lw, label=false)
end
display(plt_Vg2)
savefig(plt_Vg2,string(absPath,"/S10PazyFlutterTipMassPosRange_damp2_",idString,".pdf"))

println("Finished S10PazyFlutterTipMassPosRange.jl")