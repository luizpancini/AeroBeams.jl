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
t0, τ = 0.5, 1.0
gust = create_Continuous1DSpaceGust(spectrum=:vK,gustLength=τ*U,N=1001,σ=U/10,c0=[0, t0*U, 0],p=[0, -π/2, 0])

# Model
PazyWingContinuous1DSpaceGust,nElem,L,_ = create_Pazy(aeroSolver=aeroSolver,gustSolver=gustSolver,upright=upright,pitchAngle=θ,airspeed=U,gust=gust)

# Time variables
Δt = τ/500
tf = t0 + τ + 4t0

# Create and solve steady problem
steadyProblem = create_SteadyProblem(model=PazyWingContinuous1DSpaceGust)
solve!(steadyProblem)

# Get steady solution
x0 = steadyProblem.systemArrays.x

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=PazyWingContinuous1DSpaceGust,finalTime=tf,Δt=Δt,x0=x0)
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
relPath = "/test/outputs/figures/PazyWingContinuous1DSpaceGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=10,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingContinuous1DSpaceGust_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingContinuous1DSpaceGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingContinuous1DSpaceGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingContinuous1DSpaceGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingContinuous1DSpaceGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingContinuous1DSpaceGust_ct.pdf"))

# Normal gust velocity
W = t -> t0<t<t0+τ ? gust.W([0; t*U; 0]) : 0
plt7 = plot(xlabel="Time [s]", ylabel="Gust velocity [m/s]")
plot!(t, W.(t), color=:black, lw=lw, label=false)
display(plt7)
savefig(string(absPath,"/PazyWingContinuous1DSpaceGust_gust.pdf"))

println("Finished PazyWingContinuous1DSpaceGust.jl")

