using AeroBeams

# Flare angle [rad]
Λ = 20*π/180

# Root pitch angle [rad]
θ = 5*π/180

# Airspeed [m/s]
U = 22

# Solution method for hinge constraint
solutionMethod = "addedResidual"

# Spring stiffness [Nm/rad]
kSpring = 1e-4
kIPBendingHinge = 1e-2

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution, especially at lower airspeeds)
hasTipCorrection = true
tipLossDecayFactor = 10

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-7
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Gust
Ug = 1
ω = 5.7*2π
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2*τ,verticalVelocity=Ug)

# Model
HealySideslipFFWTOMCGustFloating = create_HealySideslipFFWT(solutionMethod=solutionMethod,flareAngle=Λ,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust) 

# Time variables
Δt = 2e-3
tf = 2

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=HealySideslipFFWTOMCGustFloating,systemSolver=NR)
solve!(steadyProblem)
println("Steady solution done")

# Create and solve dynamic problem
problem = create_DynamicProblem(model=HealySideslipFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
Nt = length(t)
tipOOP = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:Nt]
root_αₑ = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.αₑ for i in 1:Nt]
tqSpan_αₑ = [problem.aeroVariablesOverTime[i][nElementsInner].flowAnglesAndRates.αₑ for i in 1:Nt]
tip_αₑ = [problem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in 1:Nt]
root_cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:Nt]
tqSpan_cn = [problem.aeroVariablesOverTime[i][nElementsInner].aeroCoefficients.cn for i in 1:Nt]
tip_cn = [problem.aeroVariablesOverTime[i][nElem].aeroCoefficients.cn for i in 1:Nt]
root_cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:Nt]
tqSpan_cm = [problem.aeroVariablesOverTime[i][nElementsInner].aeroCoefficients.cm for i in 1:Nt]
tip_cm = [problem.aeroVariablesOverTime[i][nElem].aeroCoefficients.cm for i in 1:Nt]
root_ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:Nt]
tqSpan_ct = [problem.aeroVariablesOverTime[i][nElementsInner].aeroCoefficients.ct for i in 1:Nt]
tip_ct = [problem.aeroVariablesOverTime[i][nElem].aeroCoefficients.ct for i in 1:Nt]
tqsχ = [problem.elementalStatesOverTime[i][nElementsInner].χ for i in 1:Nt]
M2_root = [problem.nodalStatesOverTime[i][1].M_n1[2] for i in 1:Nt]

x1 = vcat([vcat(steadyProblem.model.beams[1].elements[e].x1_n1,steadyProblem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
x1_e = [steadyProblem.model.beams[1].elements[e].x1 for e in 1:nElem]
u1_e = [steadyProblem.elementalStatesOverσ[end][e].u[1] for e in 1:nElem]
p1 = vcat([vcat(steadyProblem.nodalStatesOverσ[end][e].p_n1[1],steadyProblem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
p2 = vcat([vcat(steadyProblem.nodalStatesOverσ[end][e].p_n1[2],steadyProblem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
p3 = vcat([vcat(steadyProblem.nodalStatesOverσ[end][e].p_n1[3],steadyProblem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)
steady_cn_over_span = [steadyProblem.aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:nElem]

println("Finished HealySideslipFFWTOMCGustFloating.jl")
