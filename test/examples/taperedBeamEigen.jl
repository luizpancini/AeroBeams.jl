using AeroBeams, LinearAlgebra

# Beam
L = 1
nElem = 40
EIy0 = 1
ρA0 = 1
ρA = x1 -> ρA0*(1-x1/(2*L))
EIy = x1 -> EIy0*(1-x1/(2*L))
x1 = collect(LinRange(L/nElem/2,L-L/nElem/2,nElem))
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(EIy=EIy(x)) for x in x1],I=[inertia_matrix(ρA=ρA(x)) for x in x1])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
taperedBeamEigen = create_Model(name="taperedBeamEigen",beams=[beam],BCs=[clamp])

# Create and solve the problem
nModes = 3
problem = create_EigenProblem(model=taperedBeamEigen,nModes=nModes)
solve!(problem)

# Get frequencies 
freqs = problem.frequenciesOscillatory

# Reference solution for frequencies (values of Table 3.14)
freqsRef = [4.31517; 23.5193; 63.1992]*sqrt(ρA0*L^4/EIy0)

# Show frequency comparison
ϵ_rel = freqs./freqsRef .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished taperedBeamEigen.jl")