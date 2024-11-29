using AeroBeams, LinearAlgebra

# Beam
L = 1
EIy = 1
ρA = 1
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# Spring
κ = 1
kp = [0; κ*EIy/L; 0]
spring = create_Spring(elementsIDs=[1],nodesSides=[1],kp=kp)
add_springs_to_beam!(beam=beam,springs=[spring])

# BCs
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])

# Model
pinnedSpringedBeamEigen = create_Model(name="pinnedSpringedBeamEigen",beams=[beam],BCs=[pin])

# Create and solve the problem
nModes = 2
problem = create_EigenProblem(model=pinnedSpringedBeamEigen,nModes=nModes,normalizeModeShapes=true)
solve!(problem)

# Get frequencies 
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get bending mode shapes
u3_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u3_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[3],modeShapesAbs[m].nodalStates[e].u_n2[3]) for e in 1:nElem]...)
end

# Reference solution for frequencies (values of Tables 3.10 and 3.11 for κ=1 of Hodges and Pierce (2011))
freqsRef = [1.55730; 16.2501]*sqrt(ρA*L^4/EIy)

# Show frequency comparison
ϵ_rel = freqs./freqsRef .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished pinnedSpringedBeamEigen.jl")