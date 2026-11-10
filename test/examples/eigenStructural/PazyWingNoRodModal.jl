using AeroBeams

# Meshes: :nominal and/or :refined
meshes = [:nominal, :refined]

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
problem = Array{EigenProblem}(undef,length(meshes))
eigensols = Array{EigenSolution}(undef,length(meshes))

# Loop meshes (in-vacuo wing, without the tip rod)
for (i,mesh) in enumerate(meshes)
    # Model
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),mesh=mesh,withWingtipRod=false,gravityVector=zeros(3))
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Get oscillatory frequencies in Hz
freqs = [[rads2Hz.(imag(λ)) for λ in eigensols[i].λ if imag(λ) > 0][1:nModes] for i in eachindex(eigensols)]

# Chord of the wing
_,_,chord,_ = geometrical_properties_Pazy()

# Classify the modes by their dominant displacement component, in the order of frequency (axial, in-plane bending, out-of-plane bending or torsion, the latter with the rotation weighted by half of the chord)
function mode_types(modeShapes; nModes)
    names = ["A","IP","OOP","T"]
    counts = zeros(Int,4)
    types = String[]
    for modeShape in sort(filter(m -> m.frequency > 0, modeShapes), by=m -> m.frequency)[1:nModes]
        u = reduce(hcat, [state.nodal.u2_b for state in modeShape.state])
        p = reduce(hcat, [state.nodal.p2_b for state in modeShape.state])
        energies = [sum(abs2, u[1,:]), sum(abs2, u[2,:]), sum(abs2, u[3,:]), sum(abs2, 0.5*chord*p[1,:])]
        i = argmax(energies)
        counts[i] += 1
        push!(types, string(names[i], "-", counts[i]))
    end
    return types
end
modeTypes = [mode_types(eigensols[i].modeShapes, nModes=nModes) for i in eachindex(eigensols)]

# Load reference data: natural frequencies of the Nastran shell model of the wing without the rod (SOL 103, clamped at the root, in vacuo)
modalFreqsNastran = vec(readdlm(joinpath(pkgdir(AeroBeams),"test","referenceData","PazyNoRod","S0_norod_modal_frequencies.csv"),',')[:,3])

# Nastran mode corresponding to each beam mode, matched by order: with the stringer roll tied at the ribs (nastran/tieStringerRoll.py), the shell model has no local modes, and its first modes have the same sequence of types as the beam (OOP, OOP, T, OOP, IP, OOP, T, OOP, T, OOP), so mode i of the beam is mode i of Nastran
matchedNastranFreqs = [modalFreqsNastran[1:nModes] for i in eachindex(eigensols)]

# Compute relative differences
ϵ_rel_Nastran = [freqs[i] ./ matchedNastranFreqs[i] .- 1 for i in eachindex(freqs)]

# Pretty table for modal frequency comparison
data = vcat([hcat(modeTypes[i], freqs[i], matchedNastranFreqs[i], ϵ_rel_Nastran[i] .* 100) for i in eachindex(meshes)]...)
table = pretty_table(
    String,
    data;
    column_labels    = ["Mode", "AeroBeams [Hz]", "Nastran [Hz]", "Δ Nastran [%]"],
    row_group_labels = [(i-1)*nModes + 1 => "Mesh :$(meshes[i])" for i in eachindex(meshes)],
    title            = "Modal Frequency Comparison",
    subtitle         = "Pazy wing without the tip rod: AeroBeams vs. Nastran (modes matched by order)",
    formatters       = [
        (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j == 4 && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(table)

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingNoRodModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot static mode shapes
for (i,mesh) in enumerate(meshes)
    modesPlot = plot_mode_shapes(problem[i],scale=1,nModes=nModes,modeLabels=modeTypes[i],modalColorScheme=:rainbow,view=(45,15),plotAxes=false,legendPos=(0.7,0.95),save=true,savePath=string(relPath,"/PazyWingNoRodModal_",mesh,".pdf"))
    display(modesPlot)
end

println("Finished PazyWingNoRodModal.jl")