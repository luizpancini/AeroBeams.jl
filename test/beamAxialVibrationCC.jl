using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = EA = ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=1e12,EA=EA)],I=[inertia_matrix(ρA=ρA)])

# BCs: clamped - clamped
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
beamAxialVibrationCC = create_Model(name="beamAxialVibrationCC",beams=[beam],BCs=[clamp1,clamp2])

# Create and solve the problem
nModes=6
problem = create_EigenProblem(model=beamAxialVibrationCC,nModes=nModes,getLinearSolution=true,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get torsional mode shapes
u1_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u1_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[1],modeShapesAbs[m].nodalStates[e].u_n2[1]) for e in 1:nElem]...)
end

# Analytical solution
c = sqrt(EA/ρA)
freqsAnalytical = Vector{Float64}(undef,nModes)
for m in 1:nModes
    freqsAnalytical[m] = π*c/L*m
end

# Plot
colors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_1\$")
for m in 1:nModes
    plot!(x1/L, u1_modeShapes[m], lw=2, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished beamAxialVibrationCC.jl")