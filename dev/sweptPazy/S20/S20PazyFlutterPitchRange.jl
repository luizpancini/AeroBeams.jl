using AeroBeams, DelimitedFiles, Plots, ColorSchemes

# Sweep angle [rad]
Λ = 20*π/180

# Root pitch angle range
θRange = π/180*[0,1,3,5,7,10]

# Airspeed range
URange = collect(1:1:100)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Gravity
g = 9.80665

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass configuration (assumed at tipMassPosOffset behind TE)
tipMass = 15.5e-3
tipMassPosOffset = 0.05
tipMassPos = -(chord*(1-normSparPos) + tipMassPosOffset)

# Number of modes
nModes = 5

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(θRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(θRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(θRange),nModes)
tipOOP = Array{Float64}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(θ*180/π)) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,Λ=Λ,θ=θ,airspeed=U,g=g,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0])
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

# Load reference data (experiments from Technion - https://doi.org/10.5281/zenodo.16354530)
freqs_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa0.txt")
freqs_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa1.txt")
freqs_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa3.txt")
freqs_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa5.txt")
freqs_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa7.txt")
freqs_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/aeroelastic_freqs_aoa10.txt")

# Set paths
relPath = "/dev/sweptPazy/S20/outputs/S20PazyFlutterPitchRange"
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
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,101], ylims=[0,100], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
scatter!([NaN], [NaN], ms=ms, msw=msw, mc=:white, msc=:black, label="Exp.")
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        if θ == 0 && mode <= size(freqs_aoa0_ref,1)-1
            scatter!(freqs_aoa0_ref[1,:], freqs_aoa0_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 1*π/180 && mode <= size(freqs_aoa1_ref,1)-1
            scatter!(freqs_aoa1_ref[1,:], freqs_aoa1_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 3*π/180 && mode <= size(freqs_aoa3_ref,1)-1
            scatter!(freqs_aoa3_ref[1,:], freqs_aoa3_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 5*π/180 && mode <= size(freqs_aoa5_ref,1)-1
            scatter!(freqs_aoa5_ref[1,:], freqs_aoa5_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 7*π/180 && mode <= size(freqs_aoa7_ref,1)-1
            scatter!(freqs_aoa7_ref[1,:], freqs_aoa7_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)
        elseif θ == 10*π/180 && mode <= size(freqs_aoa10_ref,1)-1
            scatter!(freqs_aoa10_ref[1,:], freqs_aoa10_ref[mode+1,:], ms=ms, msw=msw, mc=:white, msc=modeColors[mode], label=false)      
        end
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,101], ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,θ) in enumerate(θRange)
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[i,mode], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vf)
display(plt_Vg)
display(plt_Vgf)
savefig(plt_Vf,string(absPath,"/S20PazyFlutterPitchRange_freq_",tipMassConfig,".pdf"))
savefig(plt_Vg,string(absPath,"/S20PazyFlutterPitchRange_damp_",tipMassConfig,".pdf"))
savefig(plt_Vgf,string(absPath,"/S20PazyFlutterPitchRange_Vgf_",tipMassConfig,".pdf"))

println("Finished S20PazyFlutterPitchRange.jl")