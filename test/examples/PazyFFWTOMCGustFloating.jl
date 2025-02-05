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
U = 20

# Spring stiffness [Nm/rad]
kSpring = 1e-4
kIPBendingHinge = 1e1

# Solution method for hinge constraint
solutionMethod = "addedResidual"

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Gust
Ug = 1
ω = 2*2π
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=0,duration=2*τ,verticalVelocity=Ug)

# Model
PazyFFWTOMCGustFloating = create_PazyFFWT(solutionMethod=solutionMethod,hingeNode=hingeNode,flareAngle=Λ,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,foldAngle=foldAngle,gust=gust)

# Time variables
Δt = 5e-4
tf = 1.5

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=PazyFFWTOMCGustFloating,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=PazyFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x)
solve!(problem)

# Unpack numerical solution
t = problem.savedTimeVector
Nt = length(t)
tipOOP = [problem.nodalStatesOverTime[i][15].u_n2[3] for i in 1:Nt]
root_αₑ = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.αₑ for i in 1:Nt]
tqSpan_αₑ = [problem.aeroVariablesOverTime[i][hingeNode-1].flowAnglesAndRates.αₑ for i in 1:Nt]
tip_αₑ = [problem.aeroVariablesOverTime[i][15].flowAnglesAndRates.αₑ for i in 1:Nt]
root_cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:Nt]
tqSpan_cn = [problem.aeroVariablesOverTime[i][hingeNode-1].aeroCoefficients.cn for i in 1:Nt]
tip_cn = [problem.aeroVariablesOverTime[i][15].aeroCoefficients.cn for i in 1:Nt]
root_cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:Nt]
tqSpan_cm = [problem.aeroVariablesOverTime[i][hingeNode-1].aeroCoefficients.cm for i in 1:Nt]
tip_cm = [problem.aeroVariablesOverTime[i][15].aeroCoefficients.cm for i in 1:Nt]
root_ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:Nt]
tqSpan_ct = [problem.aeroVariablesOverTime[i][hingeNode-1].aeroCoefficients.ct for i in 1:Nt]
tip_ct = [problem.aeroVariablesOverTime[i][15].aeroCoefficients.ct for i in 1:Nt]
tqsχ = [problem.elementalStatesOverTime[i][hingeNode-1].χ for i in 1:Nt]
M2_root = [problem.nodalStatesOverTime[i][1].M_n1[2] for i in 1:Nt]

x1 = vcat([vcat(steadyProblem.model.beams[1].elements[e].x1_n1,steadyProblem.model.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [steadyProblem.model.beams[1].elements[e].x1 for e in 1:15]
u1_e = [steadyProblem.elementalStatesOverσ[end][e].u[1] for e in 1:15]
p1 = vcat([vcat(steadyProblem.nodalStatesOverσ[end][e].p_n1[1],steadyProblem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
p2 = vcat([vcat(steadyProblem.nodalStatesOverσ[end][e].p_n1[2],steadyProblem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
p3 = vcat([vcat(steadyProblem.nodalStatesOverσ[end][e].p_n1[3],steadyProblem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...)
steady_cn_over_span = [steadyProblem.aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]

println("Finished PazyFFWTOMCGustFloating.jl")
