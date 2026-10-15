using AeroBeams

# Root pitch angle [rad]
θ = deg2rad(5)

# Airspeed [m/s]
U = 50

# Aerodynamic and gust solvers
aeroSolver = Indicial()
gustSolver = IndicialGust(:Kussner)

# Flag for upright position
upright = true

# Gust (defined such that it begins at time t0 and lasts for τ seconds)
t0, τ = 0.1, 0.1
gust = create_DiscreteSpaceGust(type=:DARPA,gustLength=U*τ,gustWidth=2L,verticalVelocity=U*1/10,c0=[0, U*t0, 0],p=[0, -π/2, 0])

# Model
PazyWingDARPAGust,nElem,L,_ = create_Pazy(aeroSolver=aeroSolver,gustSolver=gustSolver,upright=upright,pitchAngle=θ,airspeed=U,gust=gust)

# Time variables
Δt = τ/500
tf = 10τ

# Create and solve steady problem
steadyProblem = create_SteadyProblem(model=PazyWingDARPAGust)
solve!(steadyProblem)

# Get steady solution
x0 = steadyProblem.systemArrays.x

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=PazyWingDARPAGust,finalTime=tf,Δt=Δt,x0=x0)
solve!(dynamicProblem)

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
relPath = "/test/outputs/figures/PazyWingDARPAGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=10,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingDARPAGust_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingDARPAGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingDARPAGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingDARPAGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingDARPAGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingDARPAGust_ct.pdf"))

# Gust velocity profile
s = collect(0:0.001:L)
time = collect(t0:0.01:t0+τ)
UGust(x,t) = t0<=t<=t0+τ ? (gust.UGustInertial([0;U*t;x],t))[1] : 0
X = [x for _ = time for x = s]
T = [t for t = time for _ = s]
Xnorm = X/L
Tnorm = @. (T-t0)/τ
Znorm = UGust.(X,T)/gust.verticalVelocity
plt7 = surface(ylabel="Normalized time", xlabel="Normalized span", zlabel="Normalized gust velocity")
surface!(Xnorm,Tnorm,Znorm,xlims=[0,1],ylims=[0,1],zlims=[-1,1])
display(plt7)
savefig(string(absPath,"/PazyWingDARPAGust_gust.pdf"))

println("Finished PazyWingDARPAGust.jl")