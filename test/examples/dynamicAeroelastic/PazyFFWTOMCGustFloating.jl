using AeroBeams

# Hinge node
hingeNode = 13

# Flare angle [rad]
Λ = deg2rad(10)

# Root pitch angle [rad]
θ = deg2rad(3)

# Airspeed [m/s]
U = 30

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Gust
Ug = 1
ω = Hz2rads(2)
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2τ,verticalVelocity=Ug)

# Time variables
Δt = τ/250
tf = 3

# Model
PazyFFWTOMCGustFloating,nElem,L,_ = create_Pazy(hinged=true,hingeNode=hingeNode,flareAngle=Λ,airspeed=U,pitchAngle=θ,gust=gust)

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=PazyFFWTOMCGustFloating,systemSolver=systemSolver)
solve!(steadyProblem)

# Create and solve dynamic problem
dynProblem = create_DynamicProblem(model=PazyFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x)
solve!(dynProblem)

# Unpack numerical solution
sol = dynProblem.solutionHistory
states = sol.states
hacStates = sol.hacStates
aeroVars = sol.aeroVariables
t = sol.t
ϕ = [rad2deg(hacStates[i][1].ϕ) for i in eachindex(t)]
tipOOP = [states[i][nElem].nodal.u2_b[3] for i in eachindex(t)]
M2_root = [states[i][1].nodal.M1[2] for i in eachindex(t)]
root_αₑ = [rad2deg(aeroVars[i][1].unsteadyKin.αₑ) for i in eachindex(t)]
tqSpan_αₑ = [rad2deg(aeroVars[i][round(Int,nElem*3/4)].unsteadyKin.αₑ) for i in eachindex(t)]
tip_αₑ = [rad2deg(aeroVars[i][nElem].unsteadyKin.αₑ) for i in eachindex(t)]
root_cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
tqSpan_cn = [aeroVars[i][round(Int,nElem*3/4)].aeroCoefficients.cn for i in eachindex(t)]
tip_cn = [aeroVars[i][nElem].aeroCoefficients.cn for i in eachindex(t)]
root_cm = [aeroVars[i][1].aeroCoefficients.cm for i in eachindex(t)]
tqSpan_cm = [aeroVars[i][round(Int,nElem*3/4)].aeroCoefficients.cm for i in eachindex(t)]
tip_cm = [aeroVars[i][nElem].aeroCoefficients.cm for i in eachindex(t)]
root_ct = [aeroVars[i][1].aeroCoefficients.ct for i in eachindex(t)]
tqSpan_ct = [aeroVars[i][round(Int,nElem*3/4)].aeroCoefficients.ct for i in eachindex(t)]
tip_ct = [aeroVars[i][nElem].aeroCoefficients.ct for i in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyFFWTOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, 100*tipOOP/L, c=:black, lw=lw, label=false)
display(plt_tipOOP)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_tipOOP.pdf"))

# Hinge fold angle
plt_hingeFold = plot(xlabel="Time [s]", ylabel="Hinge fold angle [deg]")
plot!(t, ϕ, c=:black, lw=lw, label=false)
display(plt_hingeFold)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_hingeFold.pdf"))

# Root OOP bending moment increment
plt_ΔM2root = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]")
plot!(t, -(M2_root .- M2_root[1]), c=:black, lw=lw, label=false)
display(plt_ΔM2root)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_ΔM2root.pdf"))

# AoA
plt_aoa = plot(xlabel="Time [s]", ylabel="Angle of attack [deg]")
plot!(t, root_αₑ, lw=lw, label="Root")
plot!(t, tqSpan_αₑ, lw=lw, label="3/4-span")
plot!(t, tip_αₑ, lw=lw, label="Tip")
display(plt_aoa)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_aoa.pdf"))

# cn
plt_cn = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, root_cn, lw=lw, label="Root")
plot!(t, tqSpan_cn, lw=lw, label="3/4-span")
plot!(t, tip_cn, lw=lw, label="Tip")
display(plt_cn)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_cn.pdf"))

# cm
plt_cm = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, root_cm, lw=lw, label="Root")
plot!(t, tqSpan_cm, lw=lw, label="3/4-span")
plot!(t, tip_cm, lw=lw, label="Tip")
display(plt_cm)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_cm.pdf"))

# ct
plt_ct = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, root_ct, lw=lw, label="Root")
plot!(t, tqSpan_ct, lw=lw, label="3/4-span")
plot!(t, tip_ct, lw=lw, label="Tip")
display(plt_ct)
savefig(string(absPath,"/PazyFFWTOMCGustFloating_ct.pdf"))

println("Finished PazyFFWTOMCGustFloating.jl")
