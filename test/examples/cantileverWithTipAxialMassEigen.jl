using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = 1
EA = 1
ρA = 1
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(EA=EA)],I=[inertia_matrix(ρA=ρA)])

# Point inertia
μ = 1
pointMass = PointInertia(elementID=nElem,mass=μ*ρA*L,η=[L/nElem/2;0;0])
add_point_inertias_to_beam!(beam,inertias=[pointMass])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
cantileverWithTipAxialMassEigen = create_Model(name="cantileverWithTipAxialMassEigen",beams=[beam],BCs=[clamp])

# Create and solve the problem
nModes = 4
problem = create_EigenProblem(model=cantileverWithTipAxialMassEigen,nModes=nModes,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs
freqsNorm = freqs*L*sqrt(ρA/EA)

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get axial mode shapes
u1_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u1_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[1],modeShapesAbs[m].nodalStates[e].u_n2[1]) for e in 1:nElem]...)
end

# Plot mode shapes
relPath = "/test/outputs/figures/cantileverWithTipAxialMassEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)
colors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_1\$")
for m in 1:nModes
    plot!(x1/L, u1_modeShapes[m], lw=2, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)
savefig(string(absPath,"/cantileverWithTipAxialMassEigen_u1.pdf"))

# Analytical solution for normalized frequencies (solutions of μ*x-cot(x) = 0 for μ=1)
freqsNormAnalytical = [0.86033; 3.4256; 6.4373; 9.5293]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished cantileverWithTipAxialMassEigen.jl")