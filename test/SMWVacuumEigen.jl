using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Wing beam
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIy,ρIz,ρIs = 0.75,5e-4,0.1-5e-4,0.1
nElem = 32
∞ = 1e12
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
SMWVacuumEigen = create_Model(name="SMWVacuumEigen",beams=[beam],BCs=[clamp])

# Create and solve problem
problem = create_EigenProblem(model=SMWVacuumEigen,nModes=5,frequencyFilterLimits=[0.01,Inf64],normalizeModeShapes=true)
solve!(problem)

# Get frequencies 
freqs = problem.frequenciesOscillatory

# Analytical solution
βL = [1.87510407; 4.69409113; 7.85475744]
FlapwiseBendingFreqsAnalytical = (βL/L).^2*sqrt(EIy/ρA)
ChordwiseBendingFreqAnalytical = (βL[1]/L).^2*sqrt(EIz/ρA)
TorsionalFreqAnalytical = π*(sqrt(GJ/(ρIs)))/L*1/2
freqsAnalytical = vcat(FlapwiseBendingFreqsAnalytical[1:2],TorsionalFreqAnalytical,ChordwiseBendingFreqAnalytical,FlapwiseBendingFreqsAnalytical[3])

# Show frequency comparison
ϵ_rel = freqs./freqsAnalytical .- 1.0
println("Relative frequency errors: $ϵ_rel")

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=5,view=(30,30),legendPos=:best,frequencyLabel="frequency",save=true,savePath="/test/outputs/figures/SMWVacuumEigen/SMWVacuumEigen_modeShapes.pdf")
display(modesPlot)

println("Finished SMWVacuumEigen.jl")