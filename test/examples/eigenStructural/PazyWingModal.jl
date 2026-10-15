using AeroBeams

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
problem = Array{EigenProblem}(undef,2)
eigensols = Array{EigenSolution}(undef,2)

# Loop position configurations
for (i,upright) in enumerate([false,true])
    # Model
    model = first(create_Pazy(upright=upright))
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,_ = mode_tracking(eigensols)

# Get frequencies of oscillatory modes in Hz
freqs,_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(), maxTracks=nModes)

# Load reference data (https://github.com/UM-A2SRL/AePW3-LDWG and https://doi.org/10.2514/1.J060621)
modalFreqsGVT = [4.39 29.80 41.00 82.50 missing; 4.26 28.50 42.00 81.50 60.70]
modalFreqsUMNAST = [4.19 28.49 41.58 83.03 100.07; 4.19 28.49 41.88 83.06 105.89]

# Compute relative differences
ϵ_rel_horz_GVT = freqs[:,1] ./ modalFreqsGVT[1,:] .- 1
ϵ_rel_vert_GVT = freqs[:,2] ./ modalFreqsGVT[2,:] .- 1
ϵ_rel_horz_UMNAST = freqs[:,1] ./ modalFreqsUMNAST[1,:] .- 1
ϵ_rel_vert_UMNAST = freqs[:,2] ./ modalFreqsUMNAST[2,:] .- 1

# Pretty table for modal frequency comparison
nModesHorz = length(ϵ_rel_horz_GVT)
nModesVert = length(ϵ_rel_vert_GVT)
data = vcat(
    hcat(["Mode $i" for i in 1:nModesHorz], freqs[1:nModesHorz, 1], modalFreqsGVT[1, :], ϵ_rel_horz_GVT .* 100,
         modalFreqsUMNAST[1, :], ϵ_rel_horz_UMNAST .* 100),
    hcat(["Mode $i" for i in 1:nModesVert]
, freqs[1:nModesVert, 2], modalFreqsGVT[2, :], ϵ_rel_vert_GVT .* 100,
         modalFreqsUMNAST[2, :], ϵ_rel_vert_UMNAST .* 100),
)
table = pretty_table(
    String,
    data;
    column_labels    = ["Mode", "AeroBeams [Hz]", "GVT [Hz]", "Δ GVT [%]", "UMNAST [Hz]", "Δ UMNAST [%]"],
    row_group_labels = [1 => "Horizontal wing", nModesHorz + 1 => "Vertical wing"],
    title            = "Modal Frequency Comparison",
    subtitle         = "AeroBeams vs. GVT and UMNAST",
    formatters       = [
        (v, i, j) -> (j in (2, 3, 5) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j in (4, 6) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(table)

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot static mode shapes
for i in [1,2]
    modesPlot = plot_mode_shapes(problem[i],scale=1,nModes=nModes,modeLabels=["OOP-1","OOP-2","T-1","OOP-3","T-IP-1"],modalColorScheme=:rainbow,view=(45,15),plotAxes=false,legendPos=(0.8,0.95),save=true,savePath=string(relPath,"/PazyWingModal_",i,".pdf"))
    display(modesPlot)
end

# Plot mode shapes animation
for i in [1,2]
    modesAnim = plot_mode_shapes_animation(problem[i],matchModeFrequency=true,nFramesPerCycle=101,nModes=nModes,modeLabels=["OOP-1","OOP-2","T-1","OOP-3","T-IP-1"],modalColorScheme=:rainbow,view=(45,15),scale=1,plotAxes=false,legendPos=(0.8,0.95),legendFontSize=10,save=true,savePath=string(relPath,"/PazyWingModal_",i,".gif"),displayProgress=true)
    display(modesAnim)
end

println("Finished PazyWingModal.jl")