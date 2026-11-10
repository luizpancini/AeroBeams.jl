using AeroBeams

# Tip mass configurations of the wing with 10 degrees of sweep
configs = [:baseline, :LE, :TE]

# Meshes: :nominal and/or :refined
meshes = [:nominal, :refined]

# Sweep angle [rad]
Λ = deg2rad(10)

# Number of elements of the wingtip rod beam (without the element of the FPS weights of the LE and TE configurations)
nElementsWingtipRod = 12

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
cases = [(config,mesh) for config in configs for mesh in meshes]
problem = Array{EigenProblem}(undef,length(cases))
eigensols = Array{EigenSolution}(undef,length(cases))
nElemWing = Vector{Int}(undef,length(cases))

# Loop configurations and meshes (in-vacuo wing, with the wingtip rod modeled as a beam)
for (i,(config,mesh)) in enumerate(cases)
    display("Solving for configuration $config, mesh $mesh")
    # Model
    model,nElemWing[i],_ = create_Pazy(aeroSolver=QuasiSteady(),sweepAngle=Λ,sweptTipMassConfig=config,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,gravityVector=zeros(3))
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Get oscillatory frequencies in Hz
freqs = [[rads2Hz.(imag(λ)) for λ in eigensols[i].λ if imag(λ) > 0][1:nModes] for i in eachindex(eigensols)]

# Chord of the wing
_,_,chord,_ = geometrical_properties_Pazy(Λ=Λ)

# Classify the modes by their dominant displacement component, in the order of frequency. Wing modes: axial (A), in-plane bending (IP), out-of-plane bending (OOP) or torsion (T), the latter with the rotation weighted by half of the chord. Rod modes (R): those in which the energy of the wing (the same displacement components, summed over its elements) is less than a fraction `rodModeRatio` of that of the rod (summed over its elements, including the FPS weights of the LE and TE configurations), as the rod follows the motion of the tip of the wing in the wing modes, but with a larger energy
rodModeRatio = 0.07
function mode_types(modeShapes; nElemWing, nModes)
    names = ["A","IP","OOP","T","R"]
    counts = zeros(Int,5)
    types = String[]
    for modeShape in sort(filter(m -> m.frequency > 0, modeShapes), by=m -> m.frequency)[1:nModes]
        u = reduce(hcat, [state.nodal.u2_b for state in modeShape.state])
        p = reduce(hcat, [state.nodal.p2_b for state in modeShape.state])
        energies(elements) = [sum(abs2, u[1,elements]), sum(abs2, u[2,elements]), sum(abs2, u[3,elements]), sum(abs2, 0.5*chord*p[1,elements])]
        wingEnergies = energies(1:nElemWing)
        rodElements = (nElemWing+2):size(u,2) # the element after those of the wing is the rigid link
        i = sum(wingEnergies) < rodModeRatio*sum(energies(rodElements)) ? 5 : argmax(wingEnergies)
        counts[i] += 1
        push!(types, string(names[i], "-", counts[i]))
    end
    return types
end
modeTypes = [mode_types(eigensols[i].modeShapes, nElemWing=nElemWing[i], nModes=nModes) for i in eachindex(eigensols)]

# Load reference data: natural frequencies of the Nastran models of the wing with the flexible tip rod (SOL 103, clamped at the root, in vacuo), with the properties of the rod and of the breaking beam of the original models, and the tip rod tied to the wing
modalFreqsNastran = [vec(readdlm(joinpath(pkgdir(AeroBeams),"test","referenceData","Pazy","S10_$(config)_modal_frequencies_flextip.csv"),',')[:,3])[1:nModes] for (config,_) in cases]

# Nastran mode corresponding to each beam mode, matched by order
matchedNastranFreqs = modalFreqsNastran

# Compute relative differences
ϵ_rel_Nastran = [freqs[i] ./ matchedNastranFreqs[i] .- 1 for i in eachindex(freqs)]

# Pretty table for modal frequency comparison
data = vcat([hcat(modeTypes[i], freqs[i], matchedNastranFreqs[i], ϵ_rel_Nastran[i] .* 100) for i in eachindex(cases)]...)
table = pretty_table(
    String,
    data;
    column_labels    = ["Mode", "AeroBeams [Hz]", "Nastran [Hz]", "Δ Nastran [%]"],
    row_group_labels = [(i-1)*nModes + 1 => "Configuration :$(cases[i][1]), mesh :$(cases[i][2])" for i in eachindex(cases)],
    title            = "Modal Frequency Comparison",
    subtitle         = "Pazy wing with 10° of sweep and the wingtip rod modeled as a beam: AeroBeams vs. Nastran (modes matched by order)",
    formatters       = [
        (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j == 4 && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(table)

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazyBeamRodModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot static mode shapes
colors = cgrad(:rainbow, nModes, categorical=true)
for (i,(config,mesh)) in enumerate(cases)
    modesPlot = plot_mode_shapes(problem[i],scale=1,nModes=nModes,modeLabels=modeTypes[i],modalColorScheme=colors,view=(45,15),plotAxes=false,legendPos=(0.7,0.95),save=true,savePath=string(relPath,"/sweptPazyBeamRodModal_",config,"_",mesh,".pdf"))
    display(modesPlot)
end

println("Finished sweptPazyBeamRodModal.jl")