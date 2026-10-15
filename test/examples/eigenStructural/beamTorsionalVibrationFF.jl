using AeroBeams

# Beam
L = GJ = ρIx = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12,GJ=GJ)],I=[inertia_matrix(ρIy=ρIx/2,ρIz=ρIx/2)])

# BCs: free - free
BCs = Vector{BC}()

# Model
beamTorsionalVibrationFF = create_Model(name="beamTorsionalVibrationFF",beams=[beam],BCs=BCs)

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes + 12

# Create and solve the problem
problem = create_EigenProblem(model=beamTorsionalVibrationFF,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter non-positive frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=1e-12, maxTracks=nModes)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Torsional mode shapes
p1_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.p1[1],
                  modeShapes[m].state[e].nodal.p2[1]] 
                  for e in 1:nElem) for m in 1:nModes]

# Analytical solution
c = sqrt(GJ/ρIx)
freqsAnalytical = [π*c/L*m for m in 1:nModes]

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/beamTorsionalVibrationFF"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
gr()

# Torsional mode shapes
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$p_1\$")
for m in 1:nModes
    plot!(x1/L, p1_modeShapes[m], lw=lw, c=colors[m], label=string("Mode $m"))
end
display(plt1)
savefig(string(absPath,"/beamTorsionalVibrationFF_p1.pdf"))

println("Finished beamTorsionalVibrationFF.jl")