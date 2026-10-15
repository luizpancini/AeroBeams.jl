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
spectrum = :vK
generationMethod = :whiteNoise
t0 = 0.5
τ = 1.0
σ = U/10
generationDuration = 200
pg = [0, -π/2, 0]
plotPSD = false
gust = create_Continuous1DGust(spectrum=spectrum,generationMethod=generationMethod,initialTime=t0,duration=τ,generationDuration=generationDuration,Uref=U,σ=σ,p=pg,plotPSD=plotPSD)

# Model
PazyWingContinuous1DGust,nElem,L,_ = create_Pazy(aeroSolver=aeroSolver,gustSolver=gustSolver,upright=upright,pitchAngle=θ,airspeed=U,gust=gust)

# Time variables
Δt = τ/500
tf = t0 + τ + 3t0

# Create and solve steady problem
steadyProblem = create_SteadyProblem(model=PazyWingContinuous1DGust)
solve!(steadyProblem)

# Get steady solution
x0 = steadyProblem.systemArrays.x

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=PazyWingContinuous1DGust,finalTime=tf,Δt=Δt,x0=x0)
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
relPath = "/test/outputs/figures/PazyWingContinuous1DGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=5,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingContinuous1DGust_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingContinuous1DGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingContinuous1DGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingContinuous1DGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingContinuous1DGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingContinuous1DGust_ct.pdf"))

# Gust velocity
W = t -> ifelse(t0<t<t0+τ,gust.W.(t),0)
plt7 = plot(xlabel="Time [s]", ylabel="Gust velocity [m/s]")
plot!(t, W.(t), color=:black, lw=lw, label=false)
display(plt7)
savefig(string(absPath,"/PazyWingContinuous1DGust_gust.pdf"))

println("Finished PazyWingContinuous1DGust.jl")

