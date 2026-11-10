using AeroBeams

# Tip mass configurations
sweptTipMassConfig = [:baseline, :LE, :TE]

# Sweep angle [rad]
Λ = deg2rad(20)

# Flag to match GVT and NASTRAN with flexible wingtip rod (true) or NASTRAN with stiffened wingtip rod (false)
flexibleWingtipRod = false

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = flexibleWingtipRod

# Flag for upright configuration
upright = true

# Gravity
applyGravity = true
g = applyGravity ? standard_atmosphere(0).g : 0.0

# Number of modes to track
nModes = 7

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 50

# Initialize outputs
problem = Array{EigenProblem}(undef,length(sweptTipMassConfig))
eigensols = Array{EigenSolution}(undef,length(sweptTipMassConfig))
freqsOsc = Array{Matrix{Float64}}(undef,length(sweptTipMassConfig))

# Loop tip mass configurations
for (i, tipMassConfig) in enumerate(sweptTipMassConfig)
    display("Solving for tipMassConfig = $tipMassConfig")
    # Model
    model,_ = create_Pazy(sweepAngle=Λ,upright=upright,gravityVector=[0,0,-g],sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)
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

# Load GVT reference data (derived from https://doi.org/10.2514/1.J066240)
GVTS20CL = [4.4, 29.6, 61.7, 86.1, missing, missing, missing]
GVTS20LE = [4.2, 28.4, 34.0, 84.4, missing, missing, missing]
GVTS20TE = [4.1, 25.1, 38.1, 83.8, missing, missing, missing]

# Load NASTRAN reference data (derived from https://doi.org/10.2514/1.J066240 with stiffened wingtip rod)
nastranS20CLF = [4.52, 29.24, 56.74, 84.65, 103.37, 161.12, 182.02] # clean, flexible (nominal) wingtip rod
nastranS20CLS = [4.52, 29.25, 57.14, 84.81, 103.44, 162.87, 188.13] # clean, stiffened wingtip rod
nastranS20LEF = [4.16, 28.16, 34.24, 83.36, 95.02, 161.08, 128.66] # LE, flexible (nominal) wingtip rod
nastranS20LES = [4.16, 28.42, 36.22, 84.14, 94.91, 166.96, 162.80] # LE, stiffened wingtip rod
nastranS20TEF = [4.13, 25.46, 38.24, 80.73, 106.44, 167.82, 186.69] # TE, flexible (nominal) wingtip rod
nastranS20TES = [4.14, 25.93, 39.86, 82.36, 94.34, 154.40, 175.77] # TE, stiffened wingtip rod

# Compute relative differences [%]
ϵ_gS20CL = 100*(freqsOsc[1] ./ GVTS20CL .- 1)
ϵ_gS20LE = 100*(freqsOsc[2] ./ GVTS20LE .- 1)
ϵ_gS20TE = 100*(freqsOsc[3] ./ GVTS20TE .- 1)
ϵ_nS20CLF = 100*(freqsOsc[1] ./ nastranS20CLF .- 1)
ϵ_nS20CLS = 100*(freqsOsc[1] ./ nastranS20CLS .- 1)
ϵ_nS20LEF = 100*(freqsOsc[2] ./ nastranS20LEF .- 1)
ϵ_nS20LES = 100*(freqsOsc[2] ./ nastranS20LES .- 1)
ϵ_nS20TEF = 100*(freqsOsc[3] ./ nastranS20TEF .- 1)
ϵ_nS20TES = 100*(freqsOsc[3] ./ nastranS20TES .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2", "OOP4"]

# Pretty tables for modal frequency comparison
dataS20CL = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[1],
    GVTS20CL,
    nastranS20CLF,
    nastranS20CLS,
    ϵ_gS20CL,
    ϵ_nS20CLF,
    ϵ_nS20CLS,
    )
dataS20LE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[2],
    GVTS20LE,
    nastranS20LEF,
    nastranS20LES,
    ϵ_gS20LE,
    ϵ_nS20LEF,
    ϵ_nS20LES,
    )
dataS20TE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[3],
    GVTS20TE,
    nastranS20TEF,
    nastranS20TES,
    ϵ_gS20TE,
    ϵ_nS20TEF,
    ϵ_nS20TES,
)

tableS20CL = pretty_table(
    String,
    dataS20CL;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S20 baseline (clean)",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20LE = pretty_table(
    String,
    dataS20LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S20 LE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20TE = pretty_table(
    String,
    dataS20TE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S20 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS20CL)
print(tableS20LE)
print(tableS20TE)

## Plots
# Set paths
relPath = "/dev/Pazy/S20/outputs/S20PazyModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
for (i,str) in enumerate(string.(sweptTipMassConfig))
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/S20PazyModal_modeShapes_$(str).pdf"))
    display(plt_swept_modeShapes)
end

println("Finished S20PazyModal.jl")