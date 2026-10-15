using AeroBeams

# Root pitch angle [rad]
θ = deg2rad(3)

# Airspeed [m/s]
U = 55

# Aerodynamic and gust solvers
aeroSolver = Indicial()
gustSolver = IndicialGust(:Kussner)

# Flag for upright position
upright = true

# Gust
σ = -1/10
ω = Hz2rads(10)
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2τ,verticalVelocity=σ*U,p=[0,-π/2,0])

# Model
PazyWingOMCGust,nElem,L,_ = create_Pazy(aeroSolver=aeroSolver,gustSolver=gustSolver,upright=upright,pitchAngle=θ,airspeed=U,gust=gust)

# System solver
# systemSolver = create_NewtonRaphson()
systemSolver = create_NLsolver(method=:newton)

# Time variables
Δt = τ/100
tf = 3

# Create and solve steady problem
steadyProblem = create_SteadyProblem(model=PazyWingOMCGust)
solve!(steadyProblem)

# Get steady solution
x0 = steadyProblem.systemArrays.x

# Create and solve dynamic problem, initialized with the steady solution
dynamicProblem = create_DynamicProblem(model=PazyWingOMCGust,systemSolver=systemSolver,finalTime=tf,Δt=Δt,x0=x0)
@time solve!(dynamicProblem)

# Unpack numerical solution
sol = dynamicProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
tipAoA = rad2deg.([aeroVars[i][nElem].unsteadyKin.αₑ for i in eachindex(t)])
tipOOP = [states[i][nElem].nodal.u2_b[3] for i in eachindex(t)]
tqSpan_cn = [aeroVars[i][12].aeroCoefficients.cn for i in eachindex(t)]
tqSpan_cm = [aeroVars[i][12].aeroCoefficients.cm for i in eachindex(t)]
tqSpan_ct = [aeroVars[i][12].aeroCoefficients.ct for i in eachindex(t)]
    
## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingOMCGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=round(Int,2e-2/Δt),plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingOMCGust_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingOMCGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingOMCGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingOMCGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingOMCGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingOMCGust_ct.pdf"))

println("Finished PazyWingOMCGust.jl")