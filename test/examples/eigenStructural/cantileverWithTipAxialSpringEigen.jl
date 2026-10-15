using AeroBeams

# Beam
L = 1
EA = 1
ρA = 1
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA)],I=[inertia_matrix(ρA=ρA)])

# Spring
μ = 1
k = μ*EA/L
spring = create_Spring(ku=k)
add_spring_to_beam!(beam=beam,nodeID=nElem+1,spring=spring)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Model
cantileverWithTipAxialSpringEigen = create_Model(name="cantileverWithTipAxialSpringEigen",beams=[beam],BCs=[clamp])

# Create and solve the problem
problem = create_EigenProblem(model=cantileverWithTipAxialSpringEigen,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter negative frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=0, maxTracks=nModes)
freqsNorm = freqs*L*sqrt(ρA/EA)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Axial mode shapes
u1_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[1],
                  modeShapes[m].state[e].nodal.u2[1]] 
                  for e in 1:nElem) for m in 1:nModes]

# Analytical solution for normalized frequencies (solutions of μ*tan(x)+x = 0 for μ=1)
freqsNormAnalytical = [2.0288; 4.9132; 7.9787; 11.086]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

# Set paths
relPath = "/test/outputs/figures/cantileverWithTipAxialSpringEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)

## Plots
# Set paths
relPath = "/test/outputs/figures/cantileverWithTipAxialSpringEigen"
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
savefig(string(absPath,"/cantileverWithTipAxialSpringEigen_u1.pdf"))

println("Finished cantileverWithTipAxialSpringEigen.jl")