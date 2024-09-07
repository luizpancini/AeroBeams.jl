using AeroBeams, ForwardDiff, DelimitedFiles

# Aerodynamic solver
aeroSolver = BLi()

# Derivation method
derivationMethod = AD()

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Frame data
airfoil = deepcopy(NACA0012)
a₀ = 0.20944
a₁ = 0.17279
b = 0.305
k = 0.098
Ma = 0.301
U = Ma*atmosphere.a

# Pitch profile
ω = k*U/b
τ = 2π/ω
t₀ = -τ/4
θ = t -> a₁*(1+sin(ω*(t+t₀)))
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Update airfoil parameters
update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=b)

# Aerodynamic surface
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=false)

# Wing
L = 1.0
EI,GJ = 2e10,2e10
nElem = 1
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=10,ρIy=1e-2,ρIz=1e-2)],rotationParametrization="E321",p0=[0;0;a₀-a₁],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0); 0.0; 0.0])

# BCs
driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t -> p(t),0,0])
journal = create_BC(name="journal",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
DSModelTest = create_Model(name="DSModelTest",beams=[wing],BCs=[driver,journal],v_A=[0;U;0])

# Set system solver options
σ0 = 1.0
maxIter = 20
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=true)

# Time variables
nCycles = 2
Δt = τ/500
tf = nCycles*τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/1e3)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=DSModelTest,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false,minΔt=Δt/2^2)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.α for i in 1:length(t)]
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:length(t)]
ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:length(t)]
cl = @. cn*cos(α) + ct*sin(α)
cdrag = @. cn*sin(α) - ct*cos(α)
Vdot3 = [problem.elementalStatesRatesOverTime[i][1].Vdot[3] for i in 1:length(t)]
χ = [problem.elementalStatesOverTime[i][1].χ for i in 1:length(t)]
χdot = [problem.elementalStatesRatesOverTime[i][1].χdot for i in 1:length(t)]

# Load reference data from McAlister et al (frame 10022)
clRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "DSModelTest", "cl.txt"))
cmRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "DSModelTest", "cm.txt"))
cdRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "DSModelTest", "cd.txt"))

println("Finished DSModelTest.jl")
