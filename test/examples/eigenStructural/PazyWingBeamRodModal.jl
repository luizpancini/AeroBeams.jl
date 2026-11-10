using AeroBeams

# Meshes: :nominal and/or :refined
meshes = [:nominal, :refined]

# Number of elements of the wingtip rod beam
nElementsWingtipRod = 16

# Flag to apply torsional correction at the root element stiffness matrix
rootTorsionCorr = true

# Number of modes to track
nModes = 15

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
models = Vector{Model}(undef,length(meshes))
problem = Array{EigenProblem}(undef,length(meshes))
eigensols = Array{EigenSolution}(undef,length(meshes))
nElemWing = Vector{Int}(undef,length(meshes))

# Loop meshes (in-vacuo wing, with the wingtip rod modeled as a beam)
for (i,mesh) in enumerate(meshes)
    # Model
    models[i],nElemWing[i],_ = create_Pazy(aeroSolver=QuasiSteady(),mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,gravityVector=zeros(3),rootTorsionCorr=rootTorsionCorr)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=models[i],eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Get oscillatory frequencies in Hz
freqs = [[rads2Hz.(imag(λ)) for λ in eigensols[i].λ if imag(λ) > 0][1:nModes] for i in eachindex(eigensols)]

# Mode shapes of the wing with the flexible tip rod
modeShapes = [eigensols[i].modeShapes[2:2:2*nModes] for i in eachindex(eigensols)]

# Chord of the wing and normalized position of the reference axis
_,_,chord,normSparPos = geometrical_properties_Pazy()

# Classify the modes by their dominant displacement component, in the order of frequency. Wing modes: axial (A), in-plane bending (IP), out-of-plane bending (OOP) or torsion (T), the latter with the rotation weighted by half of the chord. Rod modes (R): those in which the energy of the wing (the same displacement components, summed over its elements) is less than a fraction `rodModeRatio` of that of the rod (summed over its elements), as the rod follows the motion of the tip of the wing in the wing modes, but with a larger energy
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

# Load reference data: natural frequencies of the Nastran model of the wing with the flexible tip rod (SOL 103, clamped at the root, in vacuo), and its mode shapes: translations (Nastran basic frame: x aft along the chord, y along the span, z vertical) of the leading and trailing edge points of the sections of the spar, along the span from the root, and of the grids of the tip rod, along the chord from the leading edge (see extractModeShapes.py of BeamROM)
referenceDir = joinpath(pkgdir(AeroBeams),"test","referenceData","Pazy")
modalFreqsNastran = vec(readdlm(joinpath(referenceDir,"S0_modal_frequencies_flextip.csv"),',')[:,3])
wingShapesNastran = readdlm(joinpath(referenceDir,"S0_modal_shapes_flextip_wing.csv"),',',Float64,skipstart=1)
rodShapesNastran = readdlm(joinpath(referenceDir,"S0_modal_shapes_flextip_rod.csv"),',',Float64,skipstart=1)

# Linear interpolation of the data `ys` at the increasing points `xs`, evaluated at `x` (constant beyond the ends)
function interpolate_data(x,xs,ys)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs,x)
    return ys[k] + (x-xs[k])/(xs[k+1]-xs[k]) * (ys[k+1]-ys[k])
end

# Sampling points of the beam model, at the midpoints of the elements (relative to the clamp, in basis A: a1 along the span, a2 towards the leading edge and a3 vertical): those of the wing and of the rod (the element after those of the wing is the rigid link)
function sampling_elements(model,nElemWing)
    midpoints = [(model.assemblyData.r_n[n[1]] + model.assemblyData.r_n[n[2]])/2 for n in model.assemblyData.elementNodes]
    return midpoints, 1:nElemWing, (nElemWing+2):model.assemblyData.nElementsTotal
end

# Translations (in basis A) of the leading and trailing edge points of the sections of the wing and of the rod, at the sampling points, for a mode of the beam model
function beam_mode_vector(modeShape,wingElements,rodElements)
    v = Float64[]
    for e in wingElements
        u = modeShape.state[e].midpoint.u
        p = modeShape.state[e].midpoint.p
        for η in (normSparPos*chord, (normSparPos-1)*chord) # position of the leading and trailing edges along a2, from the reference axis
            append!(v, u + cross(p,[0,η,0]))
        end
    end
    for e in rodElements
        append!(v, modeShape.state[e].midpoint.u)
    end
    return v
end

# The same for a mode of the Nastran model, interpolated at the sampling points and converted to basis A (u1 = u_y, u2 = -u_x, u3 = u_z)
function nastran_mode_vector(mode,midpoints,wingElements,rodElements)
    wing = wingShapesNastran[wingShapesNastran[:,1] .== mode, :]
    rod = rodShapesNastran[rodShapesNastran[:,1] .== mode, :]
    toBasisA(u) = [u[2], -u[1], u[3]]
    v = Float64[]
    for e in wingElements, columns in (3:5, 6:8)
        append!(v, toBasisA([interpolate_data(midpoints[e][1],wing[:,2],wing[:,c]) for c in columns]))
    end
    for e in rodElements
        x = normSparPos*chord - midpoints[e][2] # chordwise position of the element, from the leading edge
        append!(v, toBasisA([interpolate_data(x,rod[:,2],rod[:,c]) for c in 3:5]))
    end
    return v
end

# Modal assurance criterion between two mode shape vectors
mac(a,b) = dot(a,b)^2 / (dot(a,a)*dot(b,b))

# Pair the modes of the beam model with those of Nastran by their shapes: the pair of largest MAC among the modes not paired yet is taken, until all the modes of the beam are paired
function pair_modes(MAC)
    pairedNastran = zeros(Int,size(MAC,1))
    availableBeam = trues(size(MAC,1))
    availableNastran = trues(size(MAC,2))
    for _ in 1:min(size(MAC)...)
        _,idx = findmax(MAC .* (availableBeam * availableNastran'))
        pairedNastran[idx[1]] = idx[2]
        availableBeam[idx[1]] = false
        availableNastran[idx[2]] = false
    end
    return pairedNastran
end

# Pairing of each mesh: Nastran mode of each beam mode, and the MAC of the pair
nastranMode = Vector{Vector{Int}}(undef,length(meshes))
macPairs = Vector{Vector{Float64}}(undef,length(meshes))
for (i,model) in enumerate(models)
    midpoints,wingElements,rodElements = sampling_elements(model,nElemWing[i])
    positiveModes = sort(filter(m -> m.frequency > 0, eigensols[i].modeShapes), by=m -> m.frequency)[1:nModes]
    beamVectors = [beam_mode_vector(m,wingElements,rodElements) for m in positiveModes]
    nastranVectors = [nastran_mode_vector(mode,midpoints,wingElements,rodElements) for mode in 1:nModes]
    MAC = [mac(a,b) for a in beamVectors, b in nastranVectors]
    nastranMode[i] = pair_modes(MAC)
    macPairs[i] = [MAC[k,nastranMode[i][k]] for k in 1:nModes]
end

# Nastran frequency corresponding to each beam mode
matchedNastranFreqs = [modalFreqsNastran[nastranMode[i]] for i in eachindex(eigensols)]

# Compute relative differences
ϵ_rel_Nastran = [freqs[i] ./ matchedNastranFreqs[i] .- 1 for i in eachindex(freqs)]

# Pretty table for modal frequency comparison
data = vcat([hcat(modeTypes[i], nastranMode[i], freqs[i], matchedNastranFreqs[i], ϵ_rel_Nastran[i] .* 100, macPairs[i]) for i in eachindex(meshes)]...)
table = pretty_table(
    String,
    data;
    column_labels    = ["Mode", "Nastran mode", "AeroBeams [Hz]", "Nastran [Hz]", "Δ Nastran [%]", "MAC"],
    row_group_labels = [(i-1)*nModes + 1 => "Mesh :$(meshes[i])" for i in eachindex(meshes)],
    title            = "Modal Frequency Comparison",
    subtitle         = "Pazy wing with the wingtip rod modeled as a beam: AeroBeams vs. Nastran (modes paired by their shapes, with the modal assurance criterion)",
    formatters       = [
        (v, i, j) -> (j in (3, 4) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j == 5 && v isa Number) ? string(round(v, digits=2), " %") : v,
        (v, i, j) -> (j == 6 && v isa Number) ? round(v, digits=3) : v,
    ],
)
print(table)

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingBeamRodModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot static mode shapes
colors = cgrad(:rainbow, nModes, categorical=true)
for (i,mesh) in enumerate(meshes)
    modesPlot = plot_mode_shapes(problem[i],scale=1,nModes=nModes,modeLabels=modeTypes[i],modalColorScheme=colors,view=(45,15),plotAxes=false,legendPos=(0.7,0.95),save=true,savePath=string(relPath,"/PazyWingBeamRodModal_",mesh,".pdf"))
    display(modesPlot)
end

# Each mode shape animation for refined mesh
for mode in 1:nModes
    modeAnim = plot_mode_shapes_animation(problem[2],modeShapes=modeShapes[2][mode:mode],modeNumbers=[mode],modalColorScheme=[colors[mode]],nFramesPerCycle=21,fps=24,scale=1,view=(45,15),legendPos=:top,plotBCs=false,plotAxes=false,displayProgress=true,save=true,savePath=string(relPath,"/PazyWingBeamRodModal_mode",mode,".gif"))
    display(modeAnim)
end

println("Finished PazyWingBeamRodModal.jl")