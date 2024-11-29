using AeroBeams

# Beam
L = EA = ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12,EA=EA)],I=[inertia_matrix(ρA=ρA)])

# BCs: free - free
BCs = Vector{BC}()

# Model
beamAxialVibrationFF = create_Model(name="beamAxialVibrationFF",beams=[beam],BCs=BCs)

# Create and solve the problem
nModes=6
problem = create_EigenProblem(model=beamAxialVibrationFF,nModes=nModes,getLinearSolution=true,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get axial mode shapes
u1_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u1_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[1],modeShapesAbs[m].nodalStates[e].u_n2[1]) for e in 1:nElem]...)
end

# Analytical solution
c = sqrt(EA/ρA)
freqsAnalytical = [π*c/L*m for m in 1:nModes]

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished beamAxialVibrationFF.jl")