using AeroBeams, LinearAlgebra

# Beam
L = 1
GJ = 1
ρIs = 1
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(GJ=GJ)],I=[inertia_matrix(ρIs=ρIs)])

# Point inertia
μ = 1
pointInertia = PointInertia(elementID=nElem,Ixx=μ*ρIs*L,η=[L/nElem/2;0;0])
add_point_inertias_to_beam!(beam,inertias=[pointInertia])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
cantileverWithTipTorsionalInertiaEigen = create_Model(name="cantileverWithTipTorsionalInertiaEigen",beams=[beam],BCs=[clamp])

# Create and solve the problem
nModes = 4
problem = create_EigenProblem(model=cantileverWithTipTorsionalInertiaEigen,nModes=nModes,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs
freqsNorm = freqs*L*sqrt(ρIs/GJ)

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get torsional mode shapes
p1_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    p1_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].p_n1[1],modeShapesAbs[m].nodalStates[e].p_n2[1]) for e in 1:nElem]...)
end

# Analytical solution for normalized frequencies (solutions of μ*x-cot(x) = 0 for μ=1)
freqsNormAnalytical = [0.86033; 3.4256; 6.4373; 9.5293]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished cantileverWithTipTorsionalInertiaEigen.jl")