using AeroBeams

# Gust frequency [Hz]
ω = 3

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Hinge configuration
hingeConfiguration = "free"

# Fold angle [rad]
foldAngle = nothing

# Flare angle [rad]
Λ = 15*π/180

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
maxIter = 200
relTol = 1e-7
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Gust
gust = create_OneMinusCosineGust(initialTime=0,duration=1/ω,verticalVelocity=Ug)

# Model
model = create_HealyBaselineFFWT(solutionMethod=solutionMethod,hingeConfiguration=hingeConfiguration,foldAngle=foldAngle,flareAngle=Λ,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,g=g,gust=gust)

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=model,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x)
solve!(problem)

# Unpack numerical solution
t = problem.savedTimeVector
M2root = [problem.nodalStatesOverTime[j][1].M_n1[2] for j in 1:length(t)]
ϕ = [problem.hingeAxisConstraintsDataOverTime[j][1].ϕ*180/π for j in 1:length(t)]

println("Finished HealyBaselineFFWTOMCGustFloating2.jl")