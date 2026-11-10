using AeroBeams

# Tip mass configurations
sweptTipMassConfig = [:baseline, :LE, :TE]

# Sweep angle [rad]
Λ = deg2rad(10)

# Mesh
mesh = :refined

# Flag to model the GVT data
matchGVT = false

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = matchGVT

# Flag for upright configuration
upright = true

# Gravity
applyGravity = matchGVT
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
    model,_ = create_Pazy(sweepAngle=Λ,mesh=mesh,upright=upright,gravityVector=[0,0,-g],sweptTipMassConfig=tipMassConfig,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)
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
GVTS10CL = [4.5, 30.2, 59.0, 87.1, missing, missing, missing]
GVTS10LE = [4.3, 28.9, 37.5, 86.2, missing, missing, 170.2]
GVTS10TE = [4.3, 27.4, 38.7, 85.2, missing, missing, 171.4]

# Load NASTRAN reference data (derived from https://doi.org/10.2514/1.J066240)
nastranS10CLF = [4.58, 29.99, 59.00, 87.01, 108.90, 167.61, 186.16] # clean, flexible wingtip rod
nastranS10CLS = [4.58, 29.99, 59.23, 87.08, 108.93, 168.48, 189.92] # clean, stiffened wingtip rod
nastranS10LEF = [4.35, 28.84, 37.83, 85.32, 103.27, 167.41, 182.49] # LE, flexible wingtip rod
nastranS10LES = [4.35, 29.12, 41.69, 86.34, 103.21, 166.67, 169.52] # LE, stiffened wingtip rod
nastranS10TEF = [4.31, 27.05, 38.71, 84.81, 104.28, 170.09, 194.99] # TE, flexible wingtip rod
nastranS10TES = [4.31, 27.64, 41.65, 85.84, 102.30, 158.63, 177.85] # TE, stiffened wingtip rod

# Compute relative differences [%]
ϵ_gS10CL = 100*(freqsOsc[1] ./ GVTS10CL .- 1)
ϵ_gS10LE = 100*(freqsOsc[2] ./ GVTS10LE .- 1)
ϵ_gS10TE = 100*(freqsOsc[3] ./ GVTS10TE .- 1)
ϵ_nS10CLS = 100*(freqsOsc[1] ./ nastranS10CLS .- 1)
ϵ_nS10CLF = 100*(freqsOsc[1] ./ nastranS10CLF .- 1)
ϵ_nS10LEF = 100*(freqsOsc[2] ./ nastranS10LEF .- 1)
ϵ_nS10LES = 100*(freqsOsc[2] ./ nastranS10LES .- 1)
ϵ_nS10TEF = 100*(freqsOsc[3] ./ nastranS10TEF .- 1)
ϵ_nS10TES = 100*(freqsOsc[3] ./ nastranS10TES .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "OOP4", "T2"]

# Pretty tables for modal frequency comparison
dataS10CL = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[1],
    GVTS10CL,
    nastranS10CLF,
    nastranS10CLS,
    ϵ_gS10CL,
    ϵ_nS10CLF,
    ϵ_nS10CLS,
    )
dataS10LE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[2],
    GVTS10LE,
    nastranS10LEF,
    nastranS10LES,
    ϵ_gS10LE,
    ϵ_nS10LEF,
    ϵ_nS10LES,
    )
dataS10TE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[3],
    GVTS10TE,
    nastranS10TEF,
    nastranS10TES,
    ϵ_gS10TE,
    ϵ_nS10TEF,
    ϵ_nS10TES,
)

tableS10CL = pretty_table(
    String,
    dataS10CL;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S10 baseline (clean)",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS10LE = pretty_table(
    String,
    dataS10LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]",
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S10 LE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS10TE = pretty_table(
    String,
    dataS10TE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]",
                        "NASTRAN (stiff) [Hz]",
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S10 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS10CL)
print(tableS10LE)
print(tableS10TE)

## Plots
# Set paths
relPath = "/dev/Pazy/S10/outputs/S10PazyModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
for (i,str) in enumerate(string.(sweptTipMassConfig))
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/S10PazyModal_modeShapes_$(str).pdf"))
    display(plt_swept_modeShapes)
end

println("Finished S10PazyModal.jl")