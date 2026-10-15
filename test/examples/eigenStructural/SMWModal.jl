using AeroBeams

# Model (in vacuo)
SMWModal,L = create_SMW(g=0,nElem=32)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Create and solve problem
problem = create_EigenProblem(model=SMWModal,eigenmethod=eigenmethod,nev=nev)
solve!(problem)

# Unpack eigenvalues
λ = problem.eigensolution.λ

# Get oscillatory frequencies and keep desired number of modes
freqs,_ = filter_tracked_modes(imag.(λ); fmin=eps(), maxTracks=nModes)

# Extract sectional properties
S = SMWModal.beams[1].sectional.S[1]
I = SMWModal.beams[1].sectional.I[1]
GJ, EIy, EIz = S[4,4], S[5,5], S[6,6]
ρA, ρIs = I[1,1], I[4,4]

# Analytical solution
βL = [1.87510407; 4.69409113; 7.85475744]

flapwiseBendingFreqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)
chordwiseBendingFreqsAnalytical = (βL/L).^2*sqrt(EIz/ρA)
torsionalFreqAnalytical = π/(2L)*sqrt(GJ/ρIs)

freqsAnalytical = vcat(flapwiseBendingFreqsAnalytical[1:2],torsionalFreqAnalytical,chordwiseBendingFreqsAnalytical[1],flapwiseBendingFreqsAnalytical[3])

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $(round.(100ϵ_rel, sigdigits=2)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,nModes=nModes,scale=5,view=(30,30),legendPos=(-0.25,0.99),save=true,savePath=string(relPath,"/SMWModal_modeShapes.pdf"))
display(modesPlot)

println("Finished SMWModal.jl")