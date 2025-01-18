using AeroBeams

# Hinge node
hingeNode = 13

# Fold angle [rad]
foldAngle = nothing

# Flare angle [rad]
Λ = 10*π/180

# Root pitch angle [rad]
θ = 3*π/180

# Airspeed [m/s]
U = 30

# Spring stiffness [Nm/rad]
kSpring = 1e-4

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-6
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Gust
Ug = -U*1/10
ω = 4*2π
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2*τ,verticalVelocity=Ug)

# Model
PazyFFWTOMCGustFloating = create_PazyFFWT(hingeNode=hingeNode,flareAngle=Λ,kSpring=kSpring,airspeed=U,pitchAngle=θ,foldAngle=foldAngle,gust=gust)

# Time variables
Δt = τ/500
tf = 70/500*τ

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=PazyFFWTOMCGustFloating,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=PazyFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x,displayFrequency=1)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tipAoA = [problem.aeroVariablesOverTime[i][15].flowAnglesAndRates.αₑ for i in 1:length(t)]
tipOOP = [problem.nodalStatesOverTime[i][15].u_n2[3] for i in 1:length(t)]
tqSpan_cn = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cn for i in 1:length(t)]
tqSpan_cm = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cm for i in 1:length(t)]
tqSpan_ct = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.ct for i in 1:length(t)]
tqsχ = [problem.elementalStatesOverTime[i][12].χ for i in 1:length(t)]

println("Finished PazyFFWTOMCGustFloating.jl")
