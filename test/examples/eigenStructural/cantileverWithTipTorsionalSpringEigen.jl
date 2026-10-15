using AeroBeams

# Beam
L = 1
GJ = 1
ρIx = 1
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(GJ=GJ)],I=[inertia_matrix(ρIy=ρIx/2,ρIz=ρIx/2)])

# Spring
μ = 1
k = μ*GJ/L
spring = create_Spring(kp=k)
add_spring_to_beam!(beam=beam,spring=spring,nodeID=nElem+1)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
cantileverWithTipTorsionalSpringEigen = create_Model(name="cantileverWithTipTorsionalSpringEigen",beams=[beam],BCs=[clamp])

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Create and solve the problem
problem = create_EigenProblem(model=cantileverWithTipTorsionalSpringEigen,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter negative frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=0, maxTracks=nModes)
freqsNorm = freqs*L*sqrt(ρIx/GJ)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Torsional mode shapes
p1_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.p1[1],
                  modeShapes[m].state[e].nodal.p2[1]] 
                  for e in 1:nElem) for m in 1:nModes]

# Analytical solution for normalized frequencies (solutions of μ*tan(x)+x = 0 for μ=1)
freqsNormAnalytical = [2.0288; 4.9132; 7.9787; 11.086]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/cantileverWithTipTorsionalSpringEigen"
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
savefig(string(absPath,"/cantileverWithTipTorsionalSpringEigen_p1.pdf"))

println("Finished cantileverWithTipTorsionalSpringEigen.jl")