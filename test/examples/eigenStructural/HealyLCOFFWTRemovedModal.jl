using AeroBeams

# Flags
wingtipRemoved = true
verticalMount = true

# Discretization
nElementsInner = 16
nElementsOuter = 6

# Number of modes to track
nModes = 9

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Model
model = create_HealyLCOFFWT(verticalMount=verticalMount,wingtipRemoved=wingtipRemoved,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter)

# Create and solve problem
problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Extract frequencies (in Hz) and mode shapes
rawFreqs = rads2Hz.(imag.(problem.eigensolution.λ))
rawDamps = real.(problem.eigensolution.λ)
rawModeShapes = problem.eigensolution.modeShapes

# Get oscillatory frequencies
freqs,_,modeShapes = filter_tracked_modes(rawFreqs,rawDamps,Vector{Union{Nothing,AeroBeams.ModeShape}}(rawModeShapes); fmin=eps(), maxTracks=nModes)

# Reference data (Fig. 5.13 of Healy's thesis)
freqsNastran = [2.6, 15.8, 15.9, 36.2, 57.2, 64.8, 99.4, 112.7, 167.8]

# Mode labels
modeLabels = ["OOP1","IP1","OOP2","OOP3","T1","OOP4","IP2","T2","OOP5"]

# Relative difference in percentage
ϵ_rel = 100 * (freqs ./ freqsNastran .- 1)

# Build table
data = hcat(modeLabels, freqs, freqsNastran, ϵ_rel)
table = pretty_table(
    String,
    data,
    column_labels = ["Mode shape", "AeroBeams [Hz]", "Nastran [Hz]", "Δ Nastran [%]"],
    formatters = [(v, i, j) -> v isa Number ? round(v, sigdigits=3) : v],
)
print(table)

# Get nodal arclength positions
nElem = nElementsInner + nElementsOuter
L = model.beams[1].geometry.ℓ
x1 = vcat([vcat(model.beams[1].elements[e].geometry.x1_n1,
                model.beams[1].elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Mode shapes
OOP_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[1],
                  modeShapes[m].state[e].nodal.u2[1]] 
                  for e in 1:nElem) for m in 1:nModes]

IP_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[2],
                  modeShapes[m].state[e].nodal.u2[2]] 
                  for e in 1:nElem) for m in 1:nModes]

T_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.p1[1],
                  modeShapes[m].state[e].nodal.p2[1]] 
                  for e in 1:nElem) for m in 1:nModes]                    

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyLCOFFWTRemovedModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
gr()

# Plot mode shapes
pltOOP = plot(ylabel="u3", title="OOP Mode Shapes")
pltIP = plot(ylabel="u2", title="IP Mode Shapes")
pltT = plot(xlabel="Nodal arclength position", ylabel="p1", title="Torsion Mode Shapes")
for m in 1:nModes
    plot!(pltOOP, x1/L, OOP_modeShapes[m], c=colors[m], lw=lw, label=false)
    plot!(pltIP, x1/L, IP_modeShapes[m], c=colors[m], lw=lw, label=false)
    plot!(pltT, x1/L, T_modeShapes[m], c=colors[m], lw=lw, label=false)
end
plt = plot(pltOOP, pltIP, pltT, layout=(3,1))
display(plt)
savefig(string(absPath,"/HealyLCOFFWTRemovedModal_modeShapes.pdf"))

println("Finished HealyLCOFFWTRemovedModal.jl")