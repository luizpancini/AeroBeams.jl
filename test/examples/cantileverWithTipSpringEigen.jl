using AeroBeams, LinearAlgebra

# Beam
L = 1
EIy = 1
ρA = 1
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# Spring
ku = [0; 0; EIy]
spring = create_Spring(elementsIDs=[nElem],nodesSides=[2],ku=ku)
add_springs_to_beam!(beam=beam,springs=[spring])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
cantileverWithTipSpringEigen = create_Model(name="cantileverWithTipSpringEigen",beams=[beam],BCs=[clamp])

# Create and solve the problem
nModes = 2
problem = create_EigenProblem(model=cantileverWithTipSpringEigen,nModes=nModes,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs
freqsNorm = freqs*L^2*sqrt(ρA/EIy)

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get bending mode shapes
u3_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u3_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[3],modeShapesAbs[m].nodalStates[e].u_n2[3]) for e in 1:nElem]...)
end

# Reference solution for normalized frequencies
freqsNormRef = [4.0399; 22.117]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormRef .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished cantileverWithTipSpringEigen.jl")