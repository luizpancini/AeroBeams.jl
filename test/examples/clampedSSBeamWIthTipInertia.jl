using AeroBeams

# Beam
L = 1
EA = 1e6
EIy = 1
ρA = 1
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# Point inertia
μ = 1
pointMass = PointInertia(elementID=nElem,mass=0,Iyy=μ*ρA*L^3,η=[L/nElem/2;0;0])
add_point_inertias_to_beam!(beam,inertias=[pointMass])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
pin = create_BC(name="pin",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])

# Model
clampedSSBeamWIthTipInertia = create_Model(name="clampedSSBeamWIthTipInertia",beams=[beam],BCs=[clamp,pin])

# Create and solve the problem
nModes = 2
problem = create_EigenProblem(model=clampedSSBeamWIthTipInertia,nModes=nModes,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get axial mode shapes
u3_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u3_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[3],modeShapesAbs[m].nodalStates[e].u_n2[3]) for e in 1:nElem]...)
end

# Reference solution for first normalized frequency at μ=1 (Problem 3.14 of Hodges and Pierce (2011))
freqRef = 1.99048*sqrt(EIy/(ρA*L^4))

# Show frequency comparison
ϵ_rel = freqs[1]/freqRef - 1.0
println("Relative error of 1st frequency: $ϵ_rel")

println("Finished clampedSSBeamWIthTipInertia.jl")