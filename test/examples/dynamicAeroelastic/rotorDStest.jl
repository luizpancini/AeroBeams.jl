using AeroBeams, ForwardDiff

# Flag to apply a transition function to airspeed and pitch profiles (suppresses noise from accelerations)
applyTransitionFunction = true

# Frame data
frame = 10_022
airfoil,a₀,a₁,b,k,Ma,U = NASA_frames_loader(frame)

# Pitch profile
ω = k*U/b
τ = 2π/ω
t₀ = -τ/4
θ = applyTransitionFunction ? 
    t -> a₁ * (1+sin(ω*(t+t₀))) * (t/τ)^10 / (1 + (t/τ)^10) : 
    t -> a₁ * (1+sin(ω*(t+t₀)))
p = t -> 4 * tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Angular velocity of the rotor (matching pitch profile)
ω3 = ω
ωRPM = rads2rpm(ω3)

# Gravity
g = standard_atmosphere(0).g

# Aerodynamic solver
aeroSolver = BLo()

# Flag to update airfoil parameters
updateAirfoilParameters = true

# Aerodynamic surface
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=updateAirfoilParameters)

# Wing
L = 10*2b
r0 = [L/4; 0; 0]
EIy,GJ = 5e6,1e8
ρA,ρIy,ρIz = 10,1e-2,1e-2
nElem = 20
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=10*EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)],p0=transform_to_WM([0,0,a₀-a₁],from=:E321),surface=surf)

# Set initial conditions
add_initial_conditions_to_beam!(beam=wing, conditions=Dict(:pdot0_of_x1 => x1->[pdot(0),0,0]))

# BCs
driver = create_BC(name="driver",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>t->p(t),:p2A=>0,:p3A=>0))

# Model
rotorDStest = create_Model(name="rotorDStest",beams=[wing],BCs=[driver],ω_A=t->[0,0,ω3],firstNodePosition=r0,gravityVector=[0;0;-g])

# System solver
NR = create_NewtonRaphson(absTol=1e-6,allowAdvanceThroughUnconvergedAeroStates=true)

# Time variables
nCycles = 4
Δt = τ/1000
tf = nCycles*τ

# Create and solve steady problem
steadyProblem = create_SteadyProblem(model=rotorDStest,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
dynProblem = create_DynamicProblem(model=rotorDStest,finalTime=tf,Δt=Δt,systemSolver=NR,x0=steadyProblem.systemArrays.x)
solve!(dynProblem)

# Unpack numerical solution
elemRangePlot = [1,15,nElem]
sol = dynProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
u3 = [[states[i][elem].nodal.u2[3] for i in eachindex(t)] for elem in elemRangePlot]
U = [[norm(aeroVars[i][elem].steadyKin.U) for i in eachindex(t)] for elem in elemRangePlot]
α = [[rad2deg.(aeroVars[i][elem].steadyKin.α) for i in eachindex(t)] for elem in elemRangePlot]
cn = [[aeroVars[i][elem].aeroCoefficients.cn for i in eachindex(t)] for elem in elemRangePlot]
cm = [[aeroVars[i][elem].aeroCoefficients.cm for i in eachindex(t)] for elem in elemRangePlot]
ct = [[aeroVars[i][elem].aeroCoefficients.ct for i in eachindex(t)] for elem in elemRangePlot]
cl = [cn[e].*cosd.(α[e]) .+ ct[e].*sind.(α[e]) for e in eachindex(elemRangePlot)]
cdrag = [cn[e].*sind.(α[e]) .- ct[e].*cosd.(α[e]) for e in eachindex(elemRangePlot)]
x1_e = [rotorDStest.assemblyData.elements[e].geometry.x1 for e in elemRangePlot]

# Load reference data from McAlister et al.
clRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cl.txt")
cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cm.txt")
cdRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cd.txt")

## Plots
# Set paths for saving figures
relPath = "/test/outputs/figures/rotorDStest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Animation
anim = plot_dynamic_deformation(dynProblem,view=(30,15),plotDistLoads=false,plotFrequency=round(Int,1e-2/Δt),plotLimits=([-L,L],[-L,L],[-L,L]),save=true,savePath=string(relPath,"/rotorDStest_deformation.gif"),displayProgress=true)
display(anim)

# Range of indices for the last cycle
rangeLastCycle = t -> (t .>= (nCycles-1)*τ) .& (t .<= nCycles*τ)

# Coefficients over time
plt_cl_t = plot(t/τ, cl, label=false, xlim=(0,nCycles), ylim=(0,2.5), ylabel="\$c_l\$")
plt_cm_t = plot(t/τ, cm, label=false, xlim=(0,nCycles), ylim=(-0.4,0.1), xlabel="Cycles", ylabel="\$c_m\$")
plt_cd_t = plot(t/τ, cdrag, label=false, xlim=(0,nCycles), ylim=(-0.1,0.7), xlabel="Cycles", ylabel="\$c_d\$")
plt1 = plot(plt_cl_t, plt_cm_t, plt_cd_t, layout=(3,1))
display(plt1)
savefig(plt1, string(absPath,"/rotorDStest_coefficients_time.png"))

# Coefficients over angle of attack
plt_cl_AoA = scatter(clRef[1,:], clRef[2,:], c=:black, ms=3, label="Exp. 2D", color=:black)
plt_cm_AoA = scatter(cmRef[1,:], cmRef[2,:], c=:black, ms=3, label=false, color=:black)
plt_cd_AoA = scatter(cdRef[1,:], cdRef[2,:], c=:black, ms=3, label=false, color=:black)
for (e,elem) in enumerate(elemRangePlot)
    plot!(plt_cl_AoA, α[e][rangeLastCycle(t)], cl[e][rangeLastCycle(t)], label="\$x/L = $(round(x1_e[e]/L,digits=2))\$", ylim=(0,2.5), xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
    plot!(plt_cm_AoA, α[e][rangeLastCycle(t)], cm[e][rangeLastCycle(t)], label=false, ylim=(-0.4,0.1), xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$")
    plot!(plt_cd_AoA, α[e][rangeLastCycle(t)], cdrag[e][rangeLastCycle(t)], label=false, ylim=(-0.1,0.7), xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
end
plt2 = plot(plt_cl_AoA, plt_cm_AoA, plt_cd_AoA, layout=(1,3))
display(plt2)
savefig(plt2, string(absPath,"/rotorDStest_coefficients_AoA.png"))

println("Finished rotorDStest.jl")