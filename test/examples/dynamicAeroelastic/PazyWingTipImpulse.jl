using AeroBeams

# Root pitch angle [rad]
θ = deg2rad(1)

# Airspeed [m/s]
U = 70

# Aerodynamic solver
aeroSolver = Indicial()

# Flag for upright position
upright = true

# Fixed geometrical properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy()

# Dummy beam for tip impulse
dummyBeam = create_Beam(length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1)])

# Set tip impulse (on dummy beam, updated later on model creation)
F₀ = 1
ω = Hz2rads(4)
τ = 2π/ω
t₀ = 0.1
F = t -> ifelse.(t.<=t₀,0.0,ifelse.(t.<=t₀+τ/2,F₀*sin.(ω*(t.-(t₀+τ))),0.0))
impulse = create_BC(name="impulse",beam=dummyBeam,node=nElem+1,conditions=Dict(:F1A=>t->F(t)))

# Model
PazyWingTipImpulse,_ = create_Pazy(aeroSolver=aeroSolver,upright=upright,pitchAngle=θ,airspeed=U,additionalBCs=[impulse])

# Time variables
Δt = τ/250
tf = 5

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=false)

# Create and solve steady problem for initial solution
steadyProblem = create_SteadyProblem(model=PazyWingTipImpulse,systemSolver=systemSolver)
solve!(steadyProblem)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=PazyWingTipImpulse,finalTime=tf,Δt=Δt,x0=steadyProblem.systemArrays.x,systemSolver=systemSolver)
@time solve!(dynamicProblem)
# @profview_allocs solve!(dynamicProblem) sample_rate=1e-3
# using Profile
# Profile.clear()
# Profile.init(n=10^7,delay=1e-4)
# @profview solve!(dynamicProblem)

# Unpack numerical solution
sol = dynamicProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
tipAoA = [rad2deg(aeroVars[i][nElem].unsteadyKin.αₑ) for i in eachindex(t)]
tipOOP = [states[i][nElem].nodal.u2_b[3] for i in eachindex(t)]
tqSpan_cn = [aeroVars[i][12].aeroCoefficients.cn for i in eachindex(t)]
tqSpan_cm = [aeroVars[i][12].aeroCoefficients.cm for i in eachindex(t)]
tqSpan_ct = [aeroVars[i][12].aeroCoefficients.ct for i in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingTipImpulse"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=round(Int,2e-2/Δt),plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingTipImpulse_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingTipImpulse_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingTipImpulse_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingTipImpulse_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingTipImpulse_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingTipImpulse_ct.pdf"))

# Root, 3/4-span and tip lift and drag coefficients over time
plot_time_outputs(dynamicProblem,elements=[1,12,nElem],elementalOutputs=["cl","cd"],save=true,saveFolder=string(relPath,"/"))

println("Finished PazyWingTipImpulse.jl")