using AeroBeams

# Hinge configuration
hingeConfiguration = "free"

# Fold angle [rad]
foldAngle = nothing

# Root pitch angle [rad]
θ = 0*π/180

# Airspeed [m/s]
U = 10

# Stiffness of the spring around the hinge for in-plane bending
kIPBendingHinge = 1e-1

# Gravity
g = 9.80665

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options
withTipCorrection = true
tipLossDecayFactor = 12

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-6
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Gust
Ug = -U*1/10
ω = 4*2π
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2*τ,verticalVelocity=Ug)

# Model
HealyBaselineFFWTOMCGustFloating = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,foldAngle=foldAngle,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,withTipCorrection=withTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,g=g,gust=gust) 

# Time variables
Δt = τ/500
tf = 5*0.2*0.07*τ

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=HealyBaselineFFWTOMCGustFloating,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=HealyBaselineFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x,displayFrequency=1)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tipAoA = [problem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in 1:length(t)]
tipOOP = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
tqSpan_cn = [problem.aeroVariablesOverTime[i][nElementsInner+1].aeroCoefficients.cn for i in 1:length(t)]
tqSpan_cm = [problem.aeroVariablesOverTime[i][nElementsInner+1].aeroCoefficients.cm for i in 1:length(t)]
tqSpan_ct = [problem.aeroVariablesOverTime[i][nElementsInner+1].aeroCoefficients.ct for i in 1:length(t)]
tqsχ = [problem.elementalStatesOverTime[i][nElementsInner+1].χ for i in 1:length(t)]

println("Finished HealyBaselineFFWTOMCGustFloating.jl")