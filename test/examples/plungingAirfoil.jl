using AeroBeams, DelimitedFiles

# Aerodynamic solver
aeroSolver = BLi()

# Atmosphere
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Test variables
airfoil = deepcopy(VERTOL23010)
k = 0.126
Ma = 0.4
θ₀ = 14.88*π/180
Δθ = 3.41*π/180
chord = 0.16
normSparPos = 1/4
U = Ma*atmosphere.a
ω = k*U/(chord/2)
hdot = t -> -U*tan(θ₀ + Δθ*sin(ω*t))

# Update airfoil parameters
update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=chord/2)

# Wing surface
surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Wing beam
L = 1
nElem = 1
∞ = 1e12
wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;0],aeroSurface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
plungingAirfoil = create_Model(name="plungingAirfoil",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere,v_A=t->[0;U;hdot(t)])

# Time variables
T = 2π/ω
cycles = 5
tf = cycles*T
Δt = T/250

# Create and solve problem
problem = create_DynamicProblem(model=plungingAirfoil,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
Nt = length(t)
α_eq = [problem.aeroVariablesOverTime[j][1].flowAnglesAndRates.α+airfoil.parametersBLi.α₀N for j in 1:Nt]
cn = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cn for j in 1:Nt]
cm = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cm for j in 1:Nt]

# Load reference data by Tyler and Leishman - Analysis of Pitch and Plunge Effects on Unsteady Airfoil Behavior (1992)
cn_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cn_exp.txt")
cn_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cn_num.txt")
cm_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cm_exp.txt")
cm_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cm_num.txt")

println("Finished plungingAirfoil.jl")