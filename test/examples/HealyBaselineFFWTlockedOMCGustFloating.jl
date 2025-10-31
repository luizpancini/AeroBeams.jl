using AeroBeams, DelimitedFiles

# Gust frequency range [Hz]
ωRange = [1,3,10]

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Hinge configuration
hingeConfiguration = "locked"

# Fold angle [rad]
foldAngle = nothing

# Root pitch angle [rad]
θ = 7.5*π/180

# Airspeed [m/s]
U = 18

# Stiffness of the spring around the hinge for in-plane bending
kIPBendingHinge = 1e4

# Gravity
g = 9.80665

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 12

# Solution method for hinge constraint
solutionMethod = "addedResidual"

# Time variables
Δt = 1e-3
tf = 2

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Initialize outputs
t = Array{Vector{Float64}}(undef,length(ωRange))
M2root = Array{Vector{Float64}}(undef,length(ωRange))
problem = Array{DynamicProblem}(undef,length(ωRange))

# Loop gust frequency range
for (i,ω) in enumerate(ωRange)
    # Display progress
    display("Solving for ω = $ω Hz")
    # Gust
    gust = create_OneMinusCosineGust(initialTime=0,duration=1/ω,verticalVelocity=Ug)
    # Model
    model = create_HealyBaselineFFWT(solutionMethod=solutionMethod,hingeConfiguration=hingeConfiguration,foldAngle=foldAngle,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,g=g,gust=gust) 
    # Solve steady problem for initial conditions
    steadyProblem = create_SteadyProblem(model=model,systemSolver=NR)
    solve!(steadyProblem)
    # Create and solve dynamic problem
    problem[i] = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x)
    solve!(problem[i])
    # Unpack numerical solution
    t[i] = problem[i].savedTimeVector
    M2root[i] = [problem[i].nodalStatesOverTime[j][1].M_n1[2] for j in 1:length(t[i])]
end

# Load reference data
ΔWRBM_Healy_1Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedOMCGustFloating/DWRBM_baseline_1Hz.txt")
ΔWRBM_Healy_3Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedOMCGustFloating/DWRBM_baseline_3Hz.txt")
ΔWRBM_Healy_10Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedOMCGustFloating/DWRBM_baseline_10Hz.txt")

println("Finished  HealyBaselineFFWTlockedOMCGustFloating.jl")