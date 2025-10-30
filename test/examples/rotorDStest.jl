using AeroBeams, ForwardDiff, DelimitedFiles

# Frame data
frame = 10_022
airfoil,a₀,a₁,b,k,Ma,U = NASA_frames_loader(frame)

# Pitch profile
ω = k*U/b
τ = 2π/ω
t₀ = -τ/4
θ = t -> a₁*(1+sin(ω*(t+t₀)))
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Angular velocity of the rotor (matching pitch profile)
ω3 = ω
ωRPM = ω3*60/(2π)

# Gravity
g = 9.80665

# Aerodynamic solver
aeroSolver = BLi()

# Flag to update airfoil parameters
updateAirfoilParameters = true

# Aerodynamic surface
surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=updateAirfoilParameters)

# Wing
L = 10*2b
r0 = [L/4; 0; 0]
EIy,GJ = 5e6,1e8
ρA,ρIy,ρIz = 10,1e-2,1e-2
nElem = 20
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=10*EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)],rotationParametrization="E321",p0=[0;0;a₀-a₁],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0); 0; 0])

# BCs
driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t->p(t),0,0])

# Model
rotorDStest = create_Model(name="rotorDStest",beams=[wing],BCs=[driver],ω_A=t->[0;0;ω3],initialPosition=r0,gravityVector=[0;0;-g])

# Set system solver options
maxIter = 50
NR = create_NewtonRaphson(maximumIterations=maxIter,displayStatus=false,allowAdvanceThroughUnconvergedAeroStates=true)

# Time variables
nCycles = 4
Δt = τ/1000
tf = nCycles*τ

# Create and solve steady problem
steadyProblem = create_SteadyProblem(model=rotorDStest,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
dynProblem = create_DynamicProblem(model=rotorDStest,finalTime=tf,Δt=Δt,systemSolver=NR,x0=steadyProblem.x)
solve!(dynProblem)

# Unpack numerical solution
elemRangePlot = [1,6,nElem]
sizeRange = 1:length(elemRangePlot)
t = dynProblem.savedTimeVector
u3 = [[dynProblem.nodalStatesOverTime[i][elem].u_n2[3] for i in 1:length(t)] for elem in elemRangePlot]
U∞rel = [[dynProblem.aeroVariablesOverTime[i][elem].flowVelocitiesAndRates.U∞ for i in 1:length(t)] for elem in elemRangePlot]
α = [[dynProblem.aeroVariablesOverTime[i][elem].flowAnglesAndRates.α for i in 1:length(t)] for elem in elemRangePlot]
cn = [[dynProblem.aeroVariablesOverTime[i][elem].aeroCoefficients.cn for i in 1:length(t)] for elem in elemRangePlot]
cm = [[dynProblem.aeroVariablesOverTime[i][elem].aeroCoefficients.cm for i in 1:length(t)] for elem in elemRangePlot]
ct = [[dynProblem.aeroVariablesOverTime[i][elem].aeroCoefficients.ct for i in 1:length(t)] for elem in elemRangePlot]
cl = @. [cn[e]*cos(α[e]) + ct[e]*sin(α[e]) for e in sizeRange]
cdrag = @. [cn[e]*sin(α[e]) - ct[e]*cos(α[e]) for e in sizeRange]

# Load reference data from McAlister et al.
clRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cl.txt")
cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cm.txt")
cdRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cd.txt")

println("Finished rotorDStest.jl")
