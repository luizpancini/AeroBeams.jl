using AeroBeams

# Beam
L = 1
EIy = 1
ρA = 1
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# Spring
μ = 1
k = μ*EIy
spring = create_Spring(ku=k)
add_spring_to_beam!(beam=beam,spring=spring,nodeID=nElem+1,p=[0;π/2;0])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
cantileverWithTipSpringEigen = create_Model(name="cantileverWithTipSpringEigen",beams=[beam],BCs=[clamp])

# System solver
solver = create_LinearSolver()

# Number of modes to track
nModes = 2

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Create and solve the problem
problem = create_EigenProblem(model=cantileverWithTipSpringEigen,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# problem = create_EigenProblem(model=cantileverWithTipSpringEigen,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
# @time solve!(problem)

# problem = create_EigenProblem(model=cantileverWithTipSpringEigen,systemSolver=solver,eigenmethod=eigenmethod,nev=nev)
# @profview_allocs solve!(problem) sample_rate=1

# using Profile
# Profile.clear()
# Profile.init(n=10^7,delay=1e-6)
# problem = create_EigenProblem(model=cantileverWithTipSpringEigen,nModes=nModes)
# @profview begin
#     for _ in 1:1000
#         solve!(problem)
#     end
# end
# @profview solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter negative frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=0, maxTracks=nModes)
freqsNorm = freqs*L^2*sqrt(ρA/EIy)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Bending mode shapes
u3_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[3],
                  modeShapes[m].state[e].nodal.u2[3]] 
                  for e in 1:nElem) for m in 1:nModes]

# Reference solution for normalized frequencies
freqsNormRef = [4.0399; 22.117]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormRef .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/cantileverWithTipSpringEigen"
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
savefig(string(absPath,"/cantileverWithTipSpringEigen_u3.pdf"))

println("Finished cantileverWithTipSpringEigen.jl")