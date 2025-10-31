using AeroBeams

# Number of modes
nModes = 5

# Gravity
g = 0

# Sweep angle [rad] range
ΛRange = vcat(0:5:30)*π/180

# Flag for ad hoc sectional stiffness corrections with sweep angle
sweepStructuralCorrections = true

# Initialize outputs
problem = Array{EigenProblem}(undef,length(ΛRange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange))

# System solver
absTol = 1e-7
NR = create_NewtonRaphson(displayStatus=true,absoluteTolerance=absTol)

# Loop sweep angle
for (i,Λ) in enumerate(ΛRange)
    display("Solving for Λ = $(round(Λ*180/π)) deg")
    # Model
    sweptPazyModal,_ = create_Pazy(Λ=Λ,upright=false,g=g,sweepStructuralCorrections=sweepStructuralCorrections)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=sweptPazyModal,nModes=nModes,systemSolver=NR)
    solve!(problem[i])
    # Get frequencies in Hz
    freqs[i] = problem[i].frequenciesOscillatory/(2π)
end

# Separate frequencies by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(ΛRange)]
end

println("Finished sweptPazyModal.jl")