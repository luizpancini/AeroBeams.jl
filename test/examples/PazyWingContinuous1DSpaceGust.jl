using AeroBeams

# Aerodynamic and gust solvers
aeroSolver = Inflow()
gustLoadsSolver = IndicialGust("Kussner")

# Flag for upright position
upright = true

# Root pitch angle [rad]
θ = 5*π/180

# Airspeed
U = 50

# Gust (defined such that it begins at time t0 and lasts for τ seconds)
spectrum = "vK"
t0 = 0.5
τ = 1.0
σ = U/10
c0 = [0;t0*U;0]
pg = [0;-π/2;0]
gust = create_Continuous1DSpaceGust(spectrum=spectrum,gustLength=τ*U,N=1001,σ=σ,c0=c0,p=pg)

# Model
PazyWingContinuous1DSpaceGust,nElem,L,_ = create_Pazy(aeroSolver=aeroSolver,gustLoadsSolver=gustLoadsSolver,upright=upright,θ=θ,airspeed=U,gust=gust)

# Set system solver options
σ0 = 1.0
maxIter = 100
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Time variables
Δt = τ/500
tf = 5*t0 + τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/10)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=PazyWingContinuous1DSpaceGust,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)

# Unpack numerical solution
t = problem.savedTimeVector
tipAoA = [problem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in 1:length(t)]
tipOOP = -[problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
tqSpan_cn = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cn for i in 1:length(t)]
tqSpan_cm = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cm for i in 1:length(t)]
tqSpan_ct = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.ct for i in 1:length(t)]
tqsχ = [problem.elementalStatesOverTime[i][12].χ for i in 1:length(t)]

println("Finished PazyWingContinuous1DSpaceGust.jl")
