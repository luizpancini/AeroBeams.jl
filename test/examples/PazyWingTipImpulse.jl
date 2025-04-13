using AeroBeams, LinearAlgebra

# Aerodynamic solver
aeroSolver = Inflow()

# Derivation method
derivationMethod = AD()

# Root pitch angle [rad]
θ = 5*π/180

# Airspeed
U = 60

# Flag for upright position
upright = true

# Fixed geometrical properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Dummy beam for tip impulse
dummyBeam = create_Beam(length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1)])

# Set tip impulse (on dummy beam, updated later on model creation)
F₀ = 10
ω = 4*2π
τ = 2π/ω
t₀ = 0.25
F = t -> ifelse.(t.<=t₀,0.0,ifelse.(t.<=t₀+τ/2,F₀*sin.(ω*(t.-(t₀+τ))),0.0))
impulse = create_BC(name="impulse",beam=dummyBeam,node=nElem+1,types=["F1A"],values=[t->F(t)])

# Model
PazyWingTipImpulse,_ = create_Pazy(aeroSolver=aeroSolver,derivationMethod=derivationMethod,upright=upright,θ=θ,airspeed=U,additionalBCs=[impulse])

# Set system solver options
σ0 = 1.0
maxIter = 100
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Time variables
Δt = τ/500
tf = t₀+5τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/10)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=PazyWingTipImpulse,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,skipInitialStatesUpdate=false)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tipAoA = [problem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in 1:length(t)]
tipOOP = -[problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
tqSpan_cn = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cn for i in 1:length(t)]
tqSpan_cm = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cm for i in 1:length(t)]
tqSpan_ct = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.ct for i in 1:length(t)]
tqsχ = [problem.elementalStatesOverTime[i][12].χ for i in 1:length(t)]
tqsχdot = [problem.elementalStatesRatesOverTime[i][12].χdot for i in 1:length(t)]

println("Finished PazyWingTipImpulse.jl")
