using AeroBeams

# Beam
L = 1
G = 1
ρ = 1
J,Is = 1,1
A,Iy,Iz = 1e6,1e6,1e6
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12,GJ=G*J)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz,ρIs=ρ*Is)])

# BCs: clamped - clamped
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
beamTorsionalVibrationCC = create_Model(name="beamTorsionalVibrationCC",beams=[beam],BCs=[clamp1,clamp2])

# Create and solve the problem
nModes=6
problem = create_EigenProblem(model=beamTorsionalVibrationCC,nModes=nModes,getLinearSolution=true,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get torsional mode shapes
p1_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    p1_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].p_n1[1],modeShapesAbs[m].nodalStates[e].p_n2[1]) for e in 1:nElem]...)
end

# Analytical solution
c = sqrt(G*J/(ρ*Is))
freqsAnalytical = [π*c/L*m for m in 1:nModes]

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished beamTorsionalVibrationCC.jl")