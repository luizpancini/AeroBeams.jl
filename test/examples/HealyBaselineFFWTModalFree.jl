using AeroBeams

# Hinge configuration
hingeConfiguration = "free"

# Discretization
nElementsInner = 32
nElementsFFWT = 8

# Gravity
g = 0

# Stiffness of the spring around the hinge for in-plane bending
kIPBendingHinge = 1e12

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-9
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Model
HealyBaselineFFWTModalFree = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,g=g,kIPBendingHinge=kIPBendingHinge,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# Number of modes
nModes = 12

# Create and solve the problem
problem = create_EigenProblem(model=HealyBaselineFFWTModalFree,nModes=nModes,frequencyFilterLimits=[1e-4,Inf],systemSolver=NR)
solve!(problem)

# Get frequencies in Hz
freqs = problem.frequenciesOscillatory/(2π)

# Reference frequencies [Hz] - see Table 3.3 of Healy's thesis
freqsRef = [0.04; 3.72; 23.24; 23.45; 65.43; 126.44; 135.01; 204.85; 257.86]

# Show frequency comparison
order = [1,2,4,5,6,8,9,10,11]
ϵ_rel = freqs[order]./freqsRef .- 1.0
println("Relative frequency errors: $ϵ_rel")

println("Finished HealyFFWTBeamModal.jl")