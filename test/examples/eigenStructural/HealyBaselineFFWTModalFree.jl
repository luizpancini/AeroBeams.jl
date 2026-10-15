using AeroBeams

# Hinge configuration
hingeConfiguration = :free

# Discretization
nElementsInner = 40
nElementsFFWT = 10

# Gravity
g = 0

# Stiffness of the spring around the hinge (used to match 1st modal frequency)
kSpring = 1.6e-4

# Model
HealyBaselineFFWTModalFree = create_HealyBaselineFFWT(aeroSolver=QuasiSteady(),hingeConfiguration=hingeConfiguration,g=g,kSpring=kSpring,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# Number of modes to track
nModes = 12

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 30

# System solver
systemSolver = create_NLsolver(jcalc=:numerical)

# Create and solve the problem
problem = create_EigenProblem(model=HealyBaselineFFWTModalFree,eigenmethod=eigenmethod,nev=nev,systemSolver=systemSolver)
solve!(problem)
# @profview_allocs solve!(problem) sample_rate=0.1
# @time solve!(problem)

# Get oscillatory frequencies
freqs,_ = filter_tracked_modes(imag.(problem.eigensolution.λ); fmin=eps(), maxTracks=nModes)

# Convert to Hz
rads2Hz!(freqs)

# Reference frequencies [Hz] - see Table 3.3 of Healy's thesis
freqsRef = [0.04, 3.72, 23.24, 23.45, 65.43, 126.44, 135.01, 204.85, 257.86]

# Compute relative differences
order = [1,2,3,4,5,7,8,10,11]
ϵ_rel = freqs[order]./freqsRef .- 1

# Build table
data = hcat(freqs[order], freqsRef, 100*ϵ_rel)
table = pretty_table(
    String,
    data,
    column_labels = ["AeroBeams [Hz]", "Healy [Hz]", "Δ [%]"],
    formatters = [(v, i, j) -> round(v, sigdigits=3)],
)
print(table)

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTModalFree"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at highest airspeed
modesPlot = plot_mode_shapes(problem,scale=1,view=(30,30),nModes=6,plotLimits=([0,1.5],[-0.75,0.75],[-0.75,0.75]),save=true,savePath=string(relPath,"/HealyFFWTBeamModal_modeShapes.pdf"))
display(modesPlot)

println("Finished HealyBaselineFFWTModalFree.jl")