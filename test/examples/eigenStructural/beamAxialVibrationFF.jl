using AeroBeams

# Beam
L = EA = ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12,EA=EA)],I=[inertia_matrix(ρA=ρA)])

# BCs: free - free
BCs = Vector{BC}()

# Model
beamAxialVibrationFF = create_Model(name="beamAxialVibrationFF",beams=[beam],BCs=BCs)

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes + 12

problem = create_EigenProblem(model=beamAxialVibrationFF,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# problem = create_EigenProblem(model=beamAxialVibrationFF,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
# @profview_allocs solve!(problem) sample_rate = 0.01

# using Profile
# Profile.clear()
# # Profile.init(n=10^7, delay=1e-4)
# problem = create_EigenProblem(model=beamAxialVibrationFF,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
# @profview solve!(problem)

# problem = create_EigenProblem(model=beamAxialVibrationFF,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
# @time solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter non-positive frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=1e-12, maxTracks=nModes)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Axial mode shapes
u1_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[1],
                  modeShapes[m].state[e].nodal.u2[1]] 
                  for e in 1:nElem) for m in 1:nModes]

# Analytical solution
c = sqrt(EA/ρA)
freqsAnalytical = [π*c/L*m for m in 1:nModes]

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/beamAxialVibrationFF"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
gr()

# Axial mode shapes
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_1\$")
for m in 1:nModes
    plot!(x1/L, u1_modeShapes[m], lw=lw, c=colors[m], label=string("Mode $m"))
end
display(plt1)
savefig(string(absPath,"/beamAxialVibrationFF_u1.pdf"))

println("Finished beamAxialVibrationFF.jl")