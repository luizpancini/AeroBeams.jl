using AeroBeams

# Beam
L = EIy = ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12,EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# BCs: free - free
BCs = Vector{BC}()

# Model
beamBendingVibrationFF = create_Model(name="beamBendingVibrationFF",beams=[beam],BCs=BCs)
# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes + 12

# Create and solve the problem
problem = create_EigenProblem(model=beamBendingVibrationFF,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter non-positive frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=1e-12, maxTracks=nModes)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Bending mode shapes
u3_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[3],
                  modeShapes[m].state[e].nodal.u2[3]] 
                  for e in 1:nElem) for m in 1:nModes]

# Analytical solution
βL = Vector{Float64}(undef,nModes)
βL[1:5] .= [4.73004074; 7.85320462; 10.9956078; 14.1371655; 17.2787597]
for m in 6:nModes
    βL[m] = (2*m+1)*π/2
end
freqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/beamBendingVibrationFF"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
gr()

# Bending mode shapes
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3\$")
for m in 1:nModes
    plot!(x1/L, u3_modeShapes[m], lw=lw, c=colors[m], label=string("Mode $m"))
end
display(plt1)
savefig(string(absPath,"/beamBendingVibrationFF_u3.pdf"))

println("Finished beamBendingVibrationFF.jl")