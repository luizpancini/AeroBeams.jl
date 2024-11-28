using AeroBeams

# Model (in vacuo)
SMWModal,L = create_SMW(g=0)

# Create and solve problem
problem = create_EigenProblem(model=SMWModal,nModes=5,frequencyFilterLimits=[0.01,Inf64],normalizeModeShapes=true)
solve!(problem)

# Get frequencies 
freqs = problem.frequenciesOscillatory

# Extract properties
GJ = SMWModal.beams[1].S[1][4,4]
EIy = SMWModal.beams[1].S[1][5,5]
EIz = SMWModal.beams[1].S[1][6,6]
ρA = SMWModal.beams[1].I[1][1,1]
ρIs = SMWModal.beams[1].I[1][4,4]

# Analytical solution
βL = [1.87510407; 4.69409113; 7.85475744]
FlapwiseBendingFreqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)
ChordwiseBendingFreqAnalytical = (βL[1]/L).^2*sqrt(EIz/ρA)
TorsionalFreqAnalytical = π*(sqrt(GJ/(ρIs)))/L*1/2
freqsAnalytical = vcat(FlapwiseBendingFreqsAnalytical[1:2],TorsionalFreqAnalytical,ChordwiseBendingFreqAnalytical,FlapwiseBendingFreqsAnalytical[3])

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished SMWModal.jl")