using AeroBeams, DelimitedFiles

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Airspeed
Ma = 0.5
U = Ma*atmosphere.a

# Wing surface data
aeroSolver = Indicial()
flapLoadsSolver = ThinAirfoilTheory()
airfoil = deepcopy(flatPlate)
chord = 0.18
normSparPos = 0.25
normFlapPos = 0.75
normFlapSpan = [0; 1]

# Flap deflection profile
k = 0.098
A = 2.5*π/180
ω = k*U/(chord/2)
δ = t -> A*sin(ω*t)

# Update airfoil parameters
update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=chord/2)

# Create wing surface
surf = create_AeroSurface(solver=aeroSolver,flapLoadsSolver=flapLoadsSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan,δ=δ)

# Wing beam
L = 1
nElem = 1
∞ = 1e10
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],aeroSurface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
flapOscillation = create_Model(name="flapOscillation",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere,v_A=[0;U;0])

# Time variables
T = 2π/ω
cycles = 10
tf = cycles*T
Δt = T/100

# Create and solve problem
problem = create_DynamicProblem(model=flapOscillation,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:length(t)]

# Load reference data by TIJDEMAN & SCHIPPERS (1973) and LEISHMAN (2006)
cnExp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "flapOscillation", "cnVsDeltaExp.txt"))
cmExp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "flapOscillation", "cmVsDeltaExp.txt"))
cnRefMod = readdlm(joinpath(dirname(@__DIR__), "referenceData", "flapOscillation", "cnVsDeltaRefMod.txt"))
cmRefMod = readdlm(joinpath(dirname(@__DIR__), "referenceData", "flapOscillation", "cmVsDeltaRefMod.txt"))

println("Finished flapOscillation.jl")