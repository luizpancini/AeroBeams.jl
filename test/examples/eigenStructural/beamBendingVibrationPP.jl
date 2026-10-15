using AeroBeams

# Beam
L = EIy = ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e6,EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# BCs: pinned - pinned
pin1 = create_BC(name="pin1",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
pin2 = create_BC(name="pin2",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))

# Model
beamBendingVibrationPP = create_Model(name="beamBendingVibrationPP",beams=[beam],BCs=[pin1,pin2])

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Create and solve the problem
problem = create_EigenProblem(model=beamBendingVibrationPP,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter negative frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=0, maxTracks=nModes)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Bending mode shapes
u3_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[3],
                  modeShapes[m].state[e].nodal.u2[3]] 
                  for e in 1:nElem) for m in 1:nModes]

# Analytical solution
βL = [m*π for m in 1:nModes]
freqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/beamBendingVibrationPP"
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
savefig(string(absPath,"/beamBendingVibrationPP_u3.pdf"))

println("Finished beamBendingVibrationPP.jl")