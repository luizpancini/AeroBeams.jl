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
L = 10*2*b
EIy,GJ = 5e6,1e8
ρA,ρIy,ρIz = 10,1e-2,1e-2
nElem = 20
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=10*EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)],rotationParametrization="E321",p0=[0;0;a₀-a₁],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0); 0.0; 0.0])

# BCs
driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t -> p(t),0,0])

# Model
wingDStest = create_Model(name="wingDStest",beams=[wing],BCs=[driver],v_A=[0;U;0])

# Set system solver options
σ0 = 1.0
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=true)

# Time variables
nCycles = 2
Δt = τ/500
tf = nCycles*τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/1e3)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=wingDStest,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false,minΔt=Δt/2^2)
solve!(problem)

# Unpack numerical solution
elemRangePlot = [1,nElem]
sizeRange = 1:length(elemRangePlot)
t = problem.timeVector
α = [[problem.aeroVariablesOverTime[i][elem].flowAnglesAndRates.α for i in 1:length(t)] for elem in elemRangePlot]
cn = [[problem.aeroVariablesOverTime[i][elem].aeroCoefficients.cn for i in 1:length(t)] for elem in elemRangePlot]
cm = [[problem.aeroVariablesOverTime[i][elem].aeroCoefficients.cm for i in 1:length(t)] for elem in elemRangePlot]
ct = [[problem.aeroVariablesOverTime[i][elem].aeroCoefficients.ct for i in 1:length(t)] for elem in elemRangePlot]
cl = @. [cn[e]*cos(α[e]) + ct[e]*sin(α[e]) for e in sizeRange]
cdrag = @. [cn[e]*sin(α[e]) - ct[e]*cos(α[e]) for e in sizeRange]

# Load reference data from McAlister et al (frame 10022)
clRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "DSModelTest", "cl.txt"))
cmRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "DSModelTest", "cm.txt"))
cdRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "DSModelTest", "cd.txt"))

println("Finished wingDStest.jl")
