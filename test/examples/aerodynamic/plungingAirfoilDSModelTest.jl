using AeroBeams

# Flag to apply a transition function to plunge profile
applyTransitionFunction = true

# Aerodynamic solver
aeroSolver = BLi()

# Atmosphere
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Test variables
airfoil = deepcopy(VERTOL23010)
k = 0.126
Ma = 0.4
θ₀ = deg2rad(14.88)
Δθ = deg2rad(3.41)
chord = 0.16
normSparPos = 1/4
U = Ma*atmosphere.a
ω = k*U/(chord/2)
T = 2π/ω
hdot = applyTransitionFunction ? 
       t -> -U * tan(θ₀ + Δθ*sin(ω*t)) * (t/T)^10 / (1 + (t/T)^10) : 
       t -> -U * tan(θ₀ + Δθ*sin(ω*t))

# Update airfoil parameters
update_all_airfoil_parameters!(airfoil,Ma)

# Wing surface
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Wing beam
L = 1
nElem = 1
wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e10)],I=[inertia_matrix(ρA=0.1)],p0=transform_to_WM([0,0,0],from=:E321),surface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))

# Model
plungingAirfoil = create_Model(name="plungingAirfoil",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere,v_A=t->[0,U,hdot(t)])

# Time variables
cycles = 10
tf = cycles*T
Δt = T/500

# System solver
systemSolver = create_NewtonRaphson(displayStatus=false)

# Create and solve problem
problem = create_DynamicProblem(model=plungingAirfoil,systemSolver=systemSolver,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
aeroVars = sol.aeroVariables
Nt = length(t)
α_eq = [rad2deg(aeroVars[i][1].steadyKin.α) for i in eachindex(t)]
αₑ = [rad2deg(aeroVars[i][1].unsteadyKin.αₑ) for i in eachindex(t)]
cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
cm = [aeroVars[i][1].aeroCoefficients.cm for i in eachindex(t)]

# Load reference data by Tyler and Leishman - Analysis of Pitch and Plunge Effects on Unsteady Airfoil Behavior (1992)
cn_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cn_exp.txt")
cn_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cn_num.txt")
cm_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cm_exp.txt")
cm_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/plungingAirfoil/cm_num.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/plungingAirfoilDSModelTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
gr()

# Index range of last cycle
rangeLastCycle = ceil(Int,(tf-T)/Δt):Nt

# cn vs α_eq
plt_cn = plot(xlabel="\$\\alpha_{eq}\$ [deg]", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs))
plot!(α_eq[rangeLastCycle], cn[rangeLastCycle], c=:black, ls=:solid, lw=lw, label=false)
plot!(cn_num[1,:], cn_num[2,:], c=:black, ls=:dashdot, lw=lw, label=false)
scatter!(cn_exp[1,:], cn_exp[2,:], c=:black, ms=ms, msw=msw, label=false)
display(plt_cn)
savefig(string(absPath,"/plungingAirfoilDSModelTest_cn.pdf"))

# cm vs α_eq
plt_cm = plot(xlabel="\$\\alpha_{eq}\$ [deg]", ylabel="\$c_m\$", tickfont=font(ts),  guidefont=font(fs), legend=:bottomleft, legendfontsize=12)
plot!(α_eq[rangeLastCycle], cm[rangeLastCycle], c=:black, ls=:solid, lw=lw, label="AeroBeams")
plot!(cm_num[1,:], cm_num[2,:], c=:black, ls=:dashdot, lw=lw, label="Model - Tyler & Leishman (1992)")
scatter!(cm_exp[1,:], cm_exp[2,:], c=:black, ms=ms, msw=msw, label="Exp. - Liiva et al. (1968)")
display(plt_cm)
savefig(string(absPath,"/plungingAirfoilDSModelTest_cm.pdf"))

println("Finished plungingAirfoilDSModelTest.jl")