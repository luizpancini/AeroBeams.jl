using AeroBeams

# Mesh: :nominal or :refined
mesh = :refined

# Flag to apply torsional correction at the root element stiffness matrix
rootTorsionCorr = true

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Upright configurations for the three tests
uprights = [false,true,true]

# Gravity for the three tests
g = [standard_atmosphere(0).g, standard_atmosphere(0).g, 0]

# Initialize outputs
problem = Array{EigenProblem}(undef,length(g))
eigensols = Array{EigenSolution}(undef,length(g))

# Loop position configurations
for (i,(upright,g)) in enumerate(zip(uprights,g))
    # Model
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),mesh=mesh,upright=upright,gravityVector=[0,0,-g],rootTorsionCorr=rootTorsionCorr)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Get oscillatory frequencies in Hz
freqs = [[rads2Hz.(imag(λ)) for λ in eigensols[i].λ if imag(λ) > 0][1:nModes] for i in eachindex(eigensols)]

# Load reference data 
modalFreqsGVT = 
[4.39 29.80 41.00 82.50 missing missing missing missing missing missing; 
 4.26 28.50 42.00 81.50 60.70 missing missing missing missing missing; 
 fill(missing, nModes)'] # (https://doi.org/10.2514/1.J060621) 
modalFreqsUMNAST = 
[4.19 28.49 41.58 83.03 100.07 missing missing missing missing missing; 
 4.19 28.49 41.88 83.06 105.89 missing missing missing missing missing; 
 fill(missing, nModes)'] # (https://github.com/UM-A2SRL/AePW3-LDWG)
modalFreqsNastran = 
[fill(missing, nModes)'; 
 fill(missing, nModes)';
 4.1934 28.5072 41.9887 83.0155 105.1676 162.6708 166.4285 278.0812 305.0567 416.8203]

# Compute relative differences
ϵ_rel_horz_GVT = freqs[1] ./ modalFreqsGVT[1,:] .- 1
ϵ_rel_vert_GVT = freqs[2] ./ modalFreqsGVT[2,:] .- 1
ϵ_rel_invacuoGVT = freqs[3] ./ modalFreqsGVT[3,:] .- 1
ϵ_rel_horz_UMNAST = freqs[1] ./ modalFreqsUMNAST[1,:] .- 1
ϵ_rel_vert_UMNAST = freqs[2] ./ modalFreqsUMNAST[2,:] .- 1
ϵ_rel_invacuo_UMNAST = freqs[3] ./ modalFreqsUMNAST[3,:] .- 1
ϵ_rel_horz_Nastran = freqs[1] ./ modalFreqsNastran[1,:] .- 1
ϵ_rel_vert_Nastran = freqs[2] ./ modalFreqsNastran[2,:] .- 1
ϵ_rel_invacuo_Nastran = freqs[3] ./ modalFreqsNastran[3,:] .- 1

# Pretty table for modal frequency comparison
data = vcat(
    hcat(["Mode $i" for i in 1:nModes], freqs[1], modalFreqsGVT[1, :], modalFreqsUMNAST[1, :], modalFreqsNastran[1, :], ϵ_rel_horz_GVT .* 100, ϵ_rel_horz_UMNAST .* 100, ϵ_rel_horz_Nastran .* 100),
    hcat(["Mode $i" for i in 1:nModes], freqs[2], modalFreqsGVT[2, :], modalFreqsUMNAST[2, :], modalFreqsNastran[2, :], ϵ_rel_vert_GVT .* 100, ϵ_rel_vert_UMNAST .* 100, ϵ_rel_vert_Nastran .* 100),
    hcat(["Mode $i" for i in 1:nModes], freqs[3], modalFreqsGVT[3, :], modalFreqsUMNAST[3, :], modalFreqsNastran[3, :], ϵ_rel_invacuoGVT .* 100, ϵ_rel_invacuo_UMNAST .* 100, ϵ_rel_invacuo_Nastran .* 100),
)
table = pretty_table(
    String,
    data;
    column_labels    = ["Mode", "AeroBeams [Hz]", "GVT [Hz]", "UM/NAST [Hz]", "Nastran [Hz]", "Δ GVT [%]", "Δ UM/NAST [%]", "Δ Nastran [%]"],
    row_group_labels = [1 => "Horizontal wing", nModes + 1 => "Vertical wing", 2*nModes + 1 => "In-vacuo wing"],
    title            = "Modal Frequency Comparison",
    subtitle         = "AeroBeams vs. GVT, UM/NAST and Nastran",
    formatters       = [
        (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(table)

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot static mode shapes
for i in [1,2,3]
    modesPlot = plot_mode_shapes(problem[i],scale=1,nModes=nModes,modeLabels=["OOP-1","OOP-2","T-1","OOP-3","T-IP-1","T-2","OOP-4","OOP-5","T-3","OOP-6"],modalColorScheme=:rainbow,view=(45,15),plotAxes=false,legendPos=(0.7,0.95),save=true,savePath=string(relPath,"/PazyWingModal_",i,".pdf"))
    display(modesPlot)
end

# Plot mode shapes animation
for i in [1,2,3]
    modesAnim = plot_mode_shapes_animation(problem[i],matchModeFrequency=true,nFramesPerCycle=101,nModes=5,modeLabels=["OOP-1","OOP-2","T-1","OOP-3","IP-1"],modalColorScheme=:rainbow,view=(45,15),scale=1,plotAxes=false,legendPos=(0.7,0.95),legendFontSize=10,save=true,savePath=string(relPath,"/PazyWingModal_",i,".gif"),displayProgress=true)
    display(modesAnim)
end

println("Finished PazyWingModal.jl")