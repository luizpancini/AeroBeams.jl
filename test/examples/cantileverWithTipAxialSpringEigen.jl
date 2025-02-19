using AeroBeams

# Beam
L = 1
EA = 1
ρA = 1
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA)],I=[inertia_matrix(ρA=ρA)])

# Spring
μ = 1
ku = [μ*EA/L; 0; 0]
spring = create_Spring(elementsIDs=[nElem],nodesSides=[2],ku=ku)
add_springs_to_beam!(beam=beam,springs=[spring])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
cantileverWithTipAxialSpringEigen = create_Model(name="cantileverWithTipAxialSpringEigen",beams=[beam],BCs=[clamp])

# Create and solve the problem
nModes = 4
problem = create_EigenProblem(model=cantileverWithTipAxialSpringEigen,nModes=nModes,normalizeModeShapes=true)
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

# Analytical solution for normalized frequencies (solutions of μ*tan(x)+x = 0 for μ=1)
freqsNormAnalytical = [2.0288; 4.9132; 7.9787; 11.086]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished cantileverWithTipAxialSpringEigen.jl")