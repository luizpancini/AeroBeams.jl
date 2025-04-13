using AeroBeams, ForwardDiff, DelimitedFiles

# Reference data source ("NASA" or "GU")
source = "NASA"

# Frame
frame = 10_022

# Number of cycles
nCycles = 4

# Starting position of the cycle ("begin" or "middle")
startingPosition = "middle"

# Flag to apply a transition function to pitch profile (suppresses noise from angular acceleration)
applyTransitionFunction = true

# Aerodynamic solver
circulatoryIndicialFunction = "Beddoes"
incompressibleInertialLoads = true
aeroSolver = BLi(circulatoryIndicialFunction=circulatoryIndicialFunction,incompressibleInertialLoads=incompressibleInertialLoads)

# Derivation method
derivationMethod = AD()

# Frame data
if source == "NASA"
    airfoil,a₀,a₁,b,k,Ma,U = NASA_frames_loader(frame)
elseif source == "GU"
    airfoil,a₀,a₁,b,k,Ma,U = GU_frames_loader(frame)
else
    error("Unavailable source")
end

# Dimensional frequency and cycle time
ω = k*U/b
τ = 2π/ω

# Pitch profile
θ = startingPosition == "begin" ? t -> a₁*(1+sin(ω*(t-τ/4))) : t -> a₁*sin(ω*t)
if applyTransitionFunction
    θ = let θ_old = θ
        t -> θ_old(t) * (t/τ)^10 / (1 + (t/τ)^10)
    end
end
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Aerodynamic surface
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=false)

# Wing (rigid)
L = 1
nElem = 1
∞ = 1e10
ρA,ρI = 10,1e-2
p0₃ = startingPosition == "begin" ? a₀-a₁ : a₀
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],rotationParametrization="E321",p0=[0;0;p0₃],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0); 0.0; 0.0])

# BCs
driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t -> p(t),0,0])
journal = create_BC(name="journal",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
pitchingAirfoilDSModelTest = create_Model(name="pitchingAirfoilDSModelTest",beams=[wing],BCs=[driver,journal],v_A=[0;U;0])

# Set system solver options
σ0 = 1.0
maxIter = 20
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=true)

# Time variables
Δt = τ/500
tf = nCycles*τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/1e2)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=pitchingAirfoilDSModelTest,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)

# Unpack numerical solution
t = problem.savedTimeVector
α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.α for i in 1:length(t)]
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:length(t)]
ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:length(t)]
cl = @. cn*cos(α) + ct*sin(α)
cdrag = @. cn*sin(α) - ct*cos(α)
χ = [problem.elementalStatesOverTime[i][1].χ for i in 1:length(t)]
χdot = [problem.elementalStatesRatesOverTime[i][1].χdot for i in 1:length(t)]

# Frame string
if source == "NASA"
    frameString = string(frame)
elseif source == "GU"
    frameString = string(frame, pad=8)
end

# Load reference data
if source == "NASA"
    clRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*frameString*"_cl.txt")
    cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*frameString*"_cm.txt")
    cdRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*frameString*"_cd.txt")
elseif source == "GU"
    cnRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/GUframes/"*frameString*"_cn.txt")
    cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/GUframes/"*frameString*"_cm.txt")
    ctRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/GUframes/"*frameString*"_ct.txt")
end

println("Finished pitchingAirfoilDSModelTest.jl")
