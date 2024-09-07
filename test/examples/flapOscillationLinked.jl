using AeroBeams

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
chord = 0.1
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

# Create wing surfaces
surf1 = create_AeroSurface(solver=aeroSolver,flapLoadsSolver=flapLoadsSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan,δ=δ)
surf2 = create_AeroSurface(solver=aeroSolver,flapLoadsSolver=flapLoadsSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan)

# Wing beams
L = 1
nElem = 1
∞ = 1e16
wing1 = create_Beam(name="wing1",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],aeroSurface=surf1)
wing2 = create_Beam(name="wing2",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],aeroSurface=surf2,connectedBeams=[wing1],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

# Flap links
m = -1
flapLink = create_FlapLink(masterBeam=wing1,slaveBeams=[wing2],δMultipliers=[m])

# BCs
clamp1 = create_BC(name="clamp1",beam=wing1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing2,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
flapOscillationLinked = create_Model(name="flapOscillationLinked",beams=[wing1,wing2],BCs=[clamp1,clamp2],flapLinks=[flapLink],atmosphere=atmosphere,v_A=[0;U;0])

# Time variables
T = 2π/ω
cycles = 5
tf = cycles*T
Δt = T/100

# Create and solve problem
problem = create_DynamicProblem(model=flapOscillationLinked,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
cnMaster = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
cmMaster = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:length(t)]
cnSlave = [problem.aeroVariablesOverTime[i][2].aeroCoefficients.cn for i in 1:length(t)]
cmSlave = [problem.aeroVariablesOverTime[i][2].aeroCoefficients.cm for i in 1:length(t)]

println("Finished flapOscillationLinked.jl")