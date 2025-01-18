using AeroBeams

# Fold angle [rad]
foldAngle = nothing

# Flare angle [rad]
Λ = 10*π/180

# Root pitch angle [rad]
θ = 3*π/180

# Airspeed [m/s]
U = 20

# Spring stiffness [Nm/rad]
kSpring = 1e-4

# Discretization
nElementsInner = 15
nElementsFFWT = 5
nElem = nElementsInner + nElementsFFWT

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution, especially at lower airspeeds)
withTipCorrection = true
tipLossDecayFactor = 10

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Gust
Ug = -U*1/10
ω = 4*2π
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2*τ,verticalVelocity=Ug)

# Model
HealyFFWTOMCGustFloating = create_HealyFFWT(flareAngle=Λ,kSpring=kSpring,airspeed=U,pitchAngle=θ,withTipCorrection=withTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust,foldAngle=foldAngle) 

# Time variables
Δt = τ/500
tf = 30/500*τ

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=HealyFFWTOMCGustFloating,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=HealyFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x,displayFrequency=1)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tipAoA = [problem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in 1:length(t)]
tipOOP = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
tqSpan_cn = [problem.aeroVariablesOverTime[i][nElementsInner+1].aeroCoefficients.cn for i in 1:length(t)]
tqSpan_cm = [problem.aeroVariablesOverTime[i][nElementsInner+1].aeroCoefficients.cm for i in 1:length(t)]
tqSpan_ct = [problem.aeroVariablesOverTime[i][nElementsInner+1].aeroCoefficients.ct for i in 1:length(t)]
tqsχ = [problem.elementalStatesOverTime[i][nElementsInner+1].χ for i in 1:length(t)]

println("Finished HealyFFWTOMCGustFloating.jl")
