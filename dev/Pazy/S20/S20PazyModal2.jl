using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(20)

# Flag to correct torsional inertia of wingtip rod (true matches GVT's, false matches NASTRAN's T1 mode) - comes from modeling as a rigid rod
sweptTipMassInertiaCorrection = true

# Tip mass configurations: baseline model + point masses removed from the LE/TE beam mass models
tipMasses, ηtipMasses, tipMassesInertias, tipMassesElementIDs = get_sweptPazy_tip_inertias_from_baseline(; Λ=Λ, configs=(:baseline, :LE, :TE), withFPSMass=true, withAdapterMass=true, withCalibrationMass=true, sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)

# Flag for upright configuration
upright = true

# Number of modes to track
nModes = 7

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 50

# Initialize outputs
problem = Array{EigenProblem}(undef,length(tipMasses))
eigensols = Array{EigenSolution}(undef,length(tipMasses))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipMasses))

# Loop tip mass configurations
for (i,m,η,I,elementIDs) in zip(1:length(tipMasses), tipMasses, ηtipMasses, tipMassesInertias, tipMassesElementIDs)
    display("Solving for configuration = $i")
    # Model
    model,_ = create_Pazy(sweepAngle=Λ,upright=upright,sweptTipMassConfig=:baseline,tipMasses=m,ηtipMasses=η,tipMassesInertias=I,tipMassesElementIDs=elementIDs)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
    # Modal frequencies
    freqs = imag.(eigensols[i].λ)
    # Oscillatory frequencies in Hz
    freqsOsc[i],_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(), maxTracks=nModes)
end

# Load GVT reference data (https://doi.org/10.2514/1.J066240)
GVTS20NO = [4.5, 29.5, 65, 86.1, missing, missing, missing]
GVTS20LE = [4.2, 28.4, 34.0, 84.4, missing, missing, missing]
GVTS20TE = [4.1, 25.1, 38.1, 83.8, missing, missing, missing]

# Load NASTRAN reference data (https://github.com/cristinariso/AePW-Swept-Pazy/tree/main/00_FEM, stiffened wingtip rod)
nastranS20NO = fill(missing, 7)
nastranS20LE = [4.16, 28.42, 36.22, 84.14, 94.91, 166.96, 162.80]
nastranS20TE = [4.14, 25.93, 39.87, 82.34, 94.34, 154.47, 175.87]

# Compute relative differences [%]
ϵ_gS20NO = 100*(freqsOsc[1] ./ GVTS20NO .- 1)
ϵ_gS20LE = 100*(freqsOsc[2] ./ GVTS20LE .- 1)
ϵ_gS20TE = 100*(freqsOsc[3] ./ GVTS20TE .- 1)
ϵ_nS20NO = 100*(freqsOsc[1] ./ nastranS20NO .- 1)
ϵ_nS20LE = 100*(freqsOsc[2] ./ nastranS20LE .- 1)
ϵ_nS20TE = 100*(freqsOsc[3] ./ nastranS20TE .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2", "OOP4"]

# Pretty tables for modal frequency comparison
dataS20NO = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[1],
    GVTS20NO,
    nastranS20NO,
    ϵ_gS20NO,
    ϵ_nS20NO,
    )
dataS20LE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[2],
    GVTS20LE,
    nastranS20LE,
    ϵ_gS20LE,
    ϵ_nS20LE,
    )
dataS20TE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[3],
    GVTS20TE,
    nastranS20TE,
    ϵ_gS20TE,
    ϵ_nS20TE,
)

tableS20NO = pretty_table(
    String,
    dataS20NO;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN [%]"],
    title            = "Modal Frequency Comparison - S20 baseline",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (5, 6) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20LE = pretty_table(
    String,
    dataS20LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN [%]"],
    title            = "Modal Frequency Comparison - S20 LE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (5, 6) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20TE = pretty_table(
    String,
    dataS20TE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN [%]"],
    title            = "Modal Frequency Comparison - S20 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (5, 6) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS20NO)
print(tableS20LE)
print(tableS20TE)

## Plots
# Set paths
relPath = "/dev/Pazy/S20/outputs/S20PazyModal2"
absPath = string(pwd(), relPath)
mkpath(absPath)

# Plot mode shapes
for i in eachindex(tipMasses)
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/S20PazyModal2_modeShapes_$(i).pdf"))
    display(plt_swept_modeShapes)
end

println("Finished S20PazyModal2.jl")