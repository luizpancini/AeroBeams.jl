using AeroBeams

# Beam
L = EIy = ρA = 1
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA)])

# Spring
μ = 1
k = μ*EIy/L
spring = create_Spring(kp=k)
add_spring_to_beam!(beam=beam,spring=spring,nodeID=1,p=[π/2;0;0])

# BCs
pin = create_BC(name="pin",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))

# Model
pinnedSpringedBeamEigen = create_Model(name="pinnedSpringedBeamEigen",beams=[beam],BCs=[pin])

# Number of modes to track
nModes = 2

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 2nModes

# Create and solve the problem
problem = create_EigenProblem(model=pinnedSpringedBeamEigen,eigenmethod=eigenmethod,nev=nev,systemSolver=systemSolver)
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

# Reference solution for frequencies (values of Tables 3.10 and 3.11 for κ=1 of Hodges and Pierce (2011))
freqsRef = [1.55730; 16.2501]*sqrt(ρA*L^4/EIy)

# Show frequency comparison
ϵ_rel = freqs./freqsRef .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/pinnedSpringedBeamEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=0.5,view=(0,0),legendPos=:outerright,save=true,savePath=string(relPath,"/pinnedSpringedBeamEigen_modeShapes.pdf"))
display(modesPlot)

println("Finished pinnedSpringedBeamEigen.jl")