using AeroBeams

# Beam
L = EIy = ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e6,EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# BCs: clamped - sliding support
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
slider = create_BC(name="slider",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
beamBendingVibrationCS = create_Model(name="beamBendingVibrationCS",beams=[beam],BCs=[clamp,slider])

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Create and solve the problem
problem = create_EigenProblem(model=beamBendingVibrationCS,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
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
βL = Vector{Float64}(undef,nModes)
βL[1:5] .= [2.36502037; 5.49780392; 8.63937983; 11.78097245; 14.92256510]
for m in 6:nModes
    βL[m] = (4*m-1)*π/4
end
freqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/beamBendingVibrationCS"
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
savefig(string(absPath,"/beamBendingVibrationCS_u3.pdf"))

println("Finished beamBendingVibrationCS.jl")