using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Aerodynamic solver
aeroSolver = BLi()

# Gust solver
gustLoadsSolver = IndicialGust("Berci&Righi")

# Derivation method
derivationMethod = AD()

# Test case (1 to 6)
testCase = 6

# Set test case data
if testCase == 1
    Ma = 0.2
    U = 68.06
    H = π
    b = 0.40663
    τ = 2*H*b/U
    w = 2.3748
    t₀ = 80*b/U
    tf = t₀ + 30*b/U
    θ = 0*π/180 
elseif testCase == 2
    Ma = 0.2
    U = 68.06
    H = π
    b = 0.40663
    τ = 2*H*b/U
    w = 2.3748
    t₀ = 80*b/U
    tf = t₀ + 30*b/U
    θ = 10*π/180  
elseif testCase == 3
    Ma = 0.2
    U = 68.06
    H = π
    b = 0.40663
    τ = 2*H*b/U
    w = 2.3748
    t₀ = 80*b/U
    tf = t₀ + 30*b/U
    θ = 15*π/180  
elseif testCase == 4
    Ma = 0.2
    U = 68.06
    H = 8π
    b = 0.40663
    τ = 2*H*b/U
    w = 2.3748
    t₀ = 80*b/U
    tf = t₀ + 80*b/U
    θ = 0*π/180
elseif testCase == 5
    Ma = 0.2
    U = 68.06
    H = 8π
    b = 0.40663
    τ = 2*H*b/U
    w = 2.3748
    t₀ = 80*b/U
    tf = t₀ + 80*b/U
    θ = 10*π/180
elseif testCase == 6
    Ma = 0.2
    U = 68.06
    H = 8π
    b = 0.40663
    τ = 2*H*b/U
    w = 2.3748
    t₀ = 60*b/U
    tf = t₀ + 80*b/U
    θ = 15*π/180         
end

# Gust
gust = create_OneMinusCosineGust(initialTime=t₀,duration=τ,verticalVelocity=w)

# Wing surface
airfoil = deepcopy(NACA0012)
surf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=2*b,normSparPos=1/4,updateAirfoilParameters=true)

# Wing beam
L = 1
nElem = 1
∞ = 1e12
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
gustTests = create_Model(name="gustTests",beams=[wing],BCs=[clamp1,clamp2],v_A=[0;U;0],gust=gust)

# Set system solver options
σ0 = 1.0
maxIter = 100
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Time variables
Δt = (tf-t₀)/1000

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=true, relaxFactor=0.5, Δt=Δt/10)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=gustTests,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false)
solve!(problem)
# @profview solve!(problem)
# @time solve!(problem)

# Unpack numerical solution
t = problem.timeVector
α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.αₑ for i in 1:length(t)]
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
fN = [problem.aeroVariablesOverTime[i][1].BLflow.fN for i in 1:length(t)]
fPrimeN = [problem.aeroVariablesOverTime[i][1].BLflow.fPrimeN for i in 1:length(t)]
f2PrimeN = [problem.elementalStatesOverTime[i][1].χ[4] for i in 1:length(t)]
r = [problem.aeroVariablesOverTime[i][1].BLkin.r for i in 1:length(t)]
cl = cn .* cos(θ)

# Load reference data
if testCase == 1
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/Hpi_A0.txt"))
elseif testCase == 2
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/Hpi_A10.txt"))
elseif testCase == 3
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/Hpi_A15.txt"))
elseif testCase == 4
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/H8pi_A0.txt"))
elseif testCase == 5
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/H8pi_A10.txt"))
elseif testCase == 6
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/H8pi_A15.txt"))
end                

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 5
i = argmin(abs.(t .- t₀))
# Separation points
plt1 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="Separation points")
plot!((t[i:end].-t[i])*U/b, fN[i:end], color=:green, lw=lw, label="\$f_N\$")
plot!((t[i:end].-t[i])*U/b, fPrimeN[i:end], color=:blue, lw=lw, label="\$f^{\\prime}_N\$")
plot!((t[i:end].-t[i])*U/b, f2PrimeN[i:end], color=:black, lw=lw, label="\$f^{\\prime\\prime}_N\$")
display(plt1)
# Δcl
plt2 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$")
plot!((t[i:end].-t[i])*U/b, cl[i:end].-cl[i], color=:black, lw=lw, label="AeroBeams")
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Mallik & Raveh (2019)")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/gustTests_dcl.pdf"))

println("Finished gustTests.jl")
