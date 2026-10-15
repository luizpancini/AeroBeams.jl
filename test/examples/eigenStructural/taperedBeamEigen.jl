using AeroBeams

# Beam
L = EIy0 = ρA0 = 1
ρA = x1 -> ρA0*(1-x1/(2*L))
EIy = x1 -> EIy0*(1-x1/(2*L))
nElem = 40
x1 = collect(LinRange(L/nElem/2,L-L/nElem/2,nElem))
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy(x)) for x in x1],I=[inertia_matrix(ρA=ρA(x)) for x in x1])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
taperedBeamEigen = create_Model(name="taperedBeamEigen",beams=[beam],BCs=[clamp])

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 30

# Create and solve the problem
problem = create_EigenProblem(model=taperedBeamEigen,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Get oscillatory frequencies and keep desired number of modes
freqs,_ = filter_tracked_modes(imag.(problem.eigensolution.λ); fmin=eps(), maxTracks=nModes)

# Reference solution for frequencies (values of Hodges & Pierce Table 3.14)
freqsRef = [4.31517; 23.5193; 63.1992]*sqrt(ρA0*L^4/EIy0)

# Show frequency comparison
ϵ_rel = freqs./freqsRef .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/taperedBeamEigen"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,nModes=nModes,scale=0.5,legendPos=:topleft,save=true,savePath=string(relPath,"/taperedBeamEigen_modeShapes.pdf"))
display(modesPlot)

println("Finished taperedBeamEigen.jl")