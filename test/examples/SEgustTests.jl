using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Aerodynamic solver
aeroSolver = BLi()

# Gust solver
gustLoadsSolver = IndicialGust("Kussner")

# Derivation method
derivationMethod = AD()

# Test case (1 to 3)
testCase = 3

# Set test case data
if testCase == 1
    Ma = 0.2
    U = 68.06
    b = 0.40663
    w = 1.1878
    t₀ = 600*b/U
    tf = t₀ + 60*b/U
    θ = 0*π/180 
elseif testCase == 2
    Ma = 0.2
    U = 68.06
    b = 0.40663
    w = 1.1878
    t₀ = 600*b/U
    tf = t₀ + 60*b/U
    θ = 10*π/180
elseif testCase == 3
    Ma = 0.2
    U = 68.06
    b = 0.40663
    w = 1.1878
    t₀ = 600*b/U
    tf = t₀ + 60*b/U
    θ = 15*π/180   
end

# Gust
gust = create_SharpEdgedGust(initialTime=t₀,verticalVelocity=w)

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
SEgustTests = create_Model(name="SEgustTests",beams=[wing],BCs=[clamp1,clamp2],v_A=[0;U;0],gust=gust)

# Set system solver options
σ0 = 1.0
maxIter = 100
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Time variables
Δt = (tf-t₀)/1000

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=true, relaxFactor=0.5, Δt=Δt/10)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=SEgustTests,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false)
solve!(problem)
# @profview solve!(problem)
# @time solve!(problem)

# Unpack numerical solution
t = problem.timeVector
α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.αₑ for i in 1:length(t)]
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:length(t)]
χ = [problem.elementalStatesOverTime[i][1].χ for i in 1:length(t)]
if typeof(aeroSolver) == BLi
    fN = [problem.aeroVariablesOverTime[i][1].BLiFlow.fN for i in 1:length(t)]
    fPrimeN = [problem.aeroVariablesOverTime[i][1].BLiFlow.fPrimeN for i in 1:length(t)]
    f2PrimeN = [problem.elementalStatesOverTime[i][1].χ[4] for i in 1:length(t)]
end
cl = @. cn*cos(θ) + ct*sin(θ)

# Load reference data
if testCase == 1
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/SE_A0.txt"))
elseif testCase == 2
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/SE_A10.txt"))
elseif testCase == 3
    clRef = readdlm(string(pwd(),"/test/referenceData/gustTests/SE_A15.txt"))
end                

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 5
i = argmin(abs.(t .- t₀))
# i=1
relPath = "/test/outputs/figures/SEgustTests"
absPath = string(pwd(),relPath)
mkpath(absPath)
gr()
# Separation points
if typeof(aeroSolver) == BLi
    plt1 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="Separation points")
    plot!((t[i:end].-t[i])*U/b, fN[i:end], color=:green, lw=lw, label="\$f_N\$")
    plot!((t[i:end].-t[i])*U/b, fPrimeN[i:end], color=:blue, lw=lw, label="\$f^{\\prime}_N\$")
    plot!((t[i:end].-t[i])*U/b, f2PrimeN[i:end], color=:black, lw=lw, label="\$f^{\\prime\\prime}_N\$")
    display(plt1)
    savefig(string(absPath,"/SEgustTests_f.pdf"))
end
# Aerodynamic states
nTotalAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nTotalAeroStates))
χ_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for j in 1:nTotalAeroStates
    χ_[j] = [χ[tt][j] for tt in 1:length(t)]
end
plt2 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="")
for j in 1:nTotalAeroStates
    plot!(t*U/b, χ_[j], c=colors[j], lw=lw, label="\$\\chi $(j)\$")
end
display(plt2)
savefig(string(absPath,"/SEgustTests_states.pdf"))
# Δcl
plt3 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$")
plot!((t[i:end].-t[i])*U/b, cl[i:end].-cl[i], color=:black, lw=lw, label="AeroBeams")
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Mallik & Raveh (2019)")
display(plt3)
savefig(string(absPath,"/SEgustTests_dcl.pdf"))

println("Finished SEgustTests.jl")
