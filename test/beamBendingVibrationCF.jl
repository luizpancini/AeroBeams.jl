using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = 1
EA,GA,GJ,EIy,EIz = 1e6,1e6,1e6,1,1e6
ρA = 1
stiffnessMatrix = diagm([EA,GA,GA,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,0,0,0])
nElem = 100
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs: clamped-free
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
beamBendingVibrationCF = create_Model(name="beamBendingVibrationCF",beams=[beam],BCs=[clamp])

# Create and solve the problem
nModes=6
problem = create_EigenProblem(model=beamBendingVibrationCF,nModes=nModes,getLinearSolution=true,normalizeModeShapes=true)
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
βL[1:5] .= [1.87510407; 4.69409113; 7.85475744; 10.99554073; 14.13716839]
for m in 6:nModes
    βL[m] = (2*m-1)*π/2
end
freqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)

# Plot
colors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
for m in 1:nModes
    plot!(x1/L, u3_modeShapes[m]/L, lw=2, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished beamBendingVibrationCF.jl")