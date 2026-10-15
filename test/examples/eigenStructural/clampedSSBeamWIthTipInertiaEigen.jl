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
pointMass = create_PointInertia(Iyy=μ*ρA*L^3,η=[L/nElem/2;0;0])
add_point_inertia_to_beam!(beam=beam,elementID=nElem,pointInertia=pointMass)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
pin = create_BC(name="pin",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))

# Model
clampedSSBeamWIthTipInertia = create_Model(name="clampedSSBeamWIthTipInertia",beams=[beam],BCs=[clamp,pin])

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Create and solve the problem
problem = create_EigenProblem(model=clampedSSBeamWIthTipInertia,eigenmethod=eigenmethod,nev=nev) 
solve!(problem)

# Unpack eigenvalues and mode shapes
λ = problem.eigensolution.λ
modeShapes = problem.eigensolution.modeShapes

# Filter negative frequencies and keep desired number of modes
freqs,_,modeShapes,_ = filter_tracked_modes(imag.(λ),real.(λ),Vector{Union{Nothing,AeroBeams.ModeShape}}(modeShapes); fmin=0, maxTracks=nModes)

# Get nodal arclength positions
x1 = vcat([vcat(beam.elements[e].geometry.x1_n1,
                beam.elements[e].geometry.x1_n2) for e in 1:nElem]...)

# Bending mode shapes
u3_modeShapes = [reduce(vcat,
                 [modeShapes[m].state[e].nodal.u1[3],
                  modeShapes[m].state[e].nodal.u2[3]] 
                  for e in 1:nElem) for m in 1:nModes]

# Reference solution for first normalized frequency at μ=1 (Problem 3.14 of Hodges and Pierce (2011))
freqRef = 1.99048*sqrt(EIy/(ρA*L^4))

# Show frequency comparison
ϵ_rel = freqs[1]/freqRef - 1.0
println("Relative error of 1st frequency: $(round(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/clampedSSBeamWIthTipInertia"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
gr()

# Plot mode shapes
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3\$")
for m in 1:nModes
    plot!(x1/L, u3_modeShapes[m], lw=lw, c=colors[m], label=string("Mode ",string(m)))
end
display(plt1)
savefig(string(absPath,"/clampedSSBeamWIthTipInertia_u3.pdf"))

println("Finished clampedSSBeamWIthTipInertiaEigen.jl")