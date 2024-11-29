using AeroBeams, LinearAlgebra

# Beam
L = 1
EA,GA,GJ,EIy,EIz = 1e6,1e6,1e6,1,1e6
ρA = 1
stiffnessMatrix = diagm([EA,GA,GA,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,0,0,0])
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix])

# BCs: clamped - clamped
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
beamBendingVibrationCC = create_Model(name="beamBendingVibrationCC",beams=[beam],BCs=[clamp1,clamp2])

# Create and solve the problem
nModes=6
problem = create_EigenProblem(model=beamBendingVibrationCC,nModes=nModes,getLinearSolution=true,normalizeModeShapes=true)
solve!(problem)

# Get frequencies and mode shapes
freqs = problem.frequenciesOscillatory
modeShapesAbs = problem.modeShapesAbs

# Get nodal arclength positions
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)

# Get bending mode shapes
u3_modeShapes = Vector{Vector{Float64}}(undef,nModes)
for m in 1:nModes
    u3_modeShapes[m] = vcat([vcat(modeShapesAbs[m].nodalStates[e].u_n1[3],modeShapesAbs[m].nodalStates[e].u_n2[3]) for e in 1:nElem]...)
end

# Analytical solution
βL = Vector{Float64}(undef,nModes)
βL[1:5] .= [4.73004074; 7.85320462; 10.9956079; 14.1371655; 17.2787597]
for m in 6:nModes
    βL[m] = (2*m+1)*π/2
end
freqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished beamBendingVibrationCC.jl")