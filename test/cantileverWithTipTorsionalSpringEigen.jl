using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = 1
GJ = 1
ρIs = 1
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(GJ=GJ)],I=[inertia_matrix(ρIs=ρIs)])

# Spring
μ = 1
kp = [μ*GJ/L; 0; 0]
spring = create_Spring(elementsIDs=[nElem],nodesSides=[2],kp=kp)
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
freqsNorm = freqs*L*sqrt(ρIs/GJ)

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get torsional mode shapes
p1_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    p1_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].p_n1[1],modeShapesAbs[m].nodalStates[e].p_n2[1]) for e in 1:nElem]...)
end

# Plot mode shapes
colors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_1\$")
for m in 1:nModes
    plot!(x1/L, p1_modeShapes[m], lw=2, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)

# Analytical solution for normalized frequencies (solutions of μ*tan(x)+x = 0 for μ=1)
freqsNormAnalytical = [2.0288; 4.9132; 7.9787; 11.086]

# Show frequency comparison
ϵ_rel = freqsNorm./freqsNormAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished cantileverWithTipAxialSpringEigen.jl")