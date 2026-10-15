using AeroBeams

# Flags
wingtipTab = true
wingtipLocked = true
verticalMount = true

# Discretization
nElementsInner = 16
nElementsOuter = 6
nElementsFFWT = 8

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Model
model = create_HealyLCOFFWT(verticalMount=verticalMount,wingtipLocked=wingtipLocked,wingtipTab=wingtipTab,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)

# Create and solve problem
problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Extract frequencies (in Hz) and mode shapes
rawFreqs = rads2Hz.(imag.(problem.eigensolution.λ))
rawDamps = real.(problem.eigensolution.λ)
rawModeShapes = problem.eigensolution.modeShapes

# Get oscillatory frequencies
freqs,_,modeShapes = filter_tracked_modes(rawFreqs,rawDamps,Vector{Union{Nothing,AeroBeams.ModeShape}}(rawModeShapes); fmin=eps(), maxTracks=nModes)

# Reference data (Tables 5.3 and 5.4 of Healy's thesis)
freqsGVT = [1.62, 8.81, 9.31, 22.37, 30.53, 30.74]
freqsNastran = freqsGVT .* [1.03, 0.9, 0.8, 0.92, 0.97, 0.63]

# Mode labels
modeLabels = ["OOP1","OOP2","IP1","OOP3","IP2","T"]

# Relative errors
ϵ_rel_GVT = 100 * (freqs ./ freqsGVT .- 1)
ϵ_rel_Nastran = 100 * (freqs ./ freqsNastran .- 1)

# Build table
data = hcat(modeLabels, freqs, freqsGVT, freqsNastran, ϵ_rel_GVT, ϵ_rel_Nastran)
table = pretty_table(
    String,
    data,
    column_labels = ["Mode shape", "AeroBeams [Hz]", "GVT [Hz]", "Nastran [Hz]", "Δ GVT [%]", "Δ Nastran [%]"],
    formatters = [(v, i, j) -> v isa Number ? round(v, sigdigits=3) : v],
)
print(table)

# Get nodal arclength positions
nElem = nElementsInner + nElementsOuter + nElementsFFWT
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
relPath = "/test/outputs/figures/HealyLCOFFWTLockedModal"
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
savefig(string(absPath,"/HealyLCOFFWTLockedModal_modeShapes.pdf"))

println("Finished HealyLCOFFWTLockedModal.jl")