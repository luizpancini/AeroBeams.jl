using AeroBeams

# Flare angle [rad]
Λ = deg2rad(20)

# Root pitch angle [rad]
θ = deg2rad(5)

# Airspeed [m/s]
U = 22

# Discretization
nElementsInner = 15
nElementsFFWT = 5
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 10

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Gust
Ug = 1
ω = Hz2rads(5.7)
τ = 2π/ω
gust = create_OneMinusCosineGust(initialTime=τ,duration=2τ,verticalVelocity=Ug)

# Model
HealySideslipFFWTOMCGustFloating = create_HealySideslipFFWT(flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust)

# Fixed geometry
beam = HealySideslipFFWTOMCGustFloating.beams[1]
L = beam.geometry.ℓ
x1_n = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                    beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)
x1_e = [beam.elements[e].geometry.x1 for e in 1:nElem]

# Time variables
Δt = 2.5e-3
tf = 2

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=HealySideslipFFWTOMCGustFloating,systemSolver=systemSolver)
solve!(steadyProblem)

# Steady solution outputs
steadySol = steadyProblem.loadstepHistory
steadyState = steadySol.states[end]
steadyAeroVars = steadySol.aeroVariables[end]
p1s = reduce(vcat,[steadyState[e].nodal.p1[1],
                   steadyState[e].nodal.p2[1]] for e in 1:nElem)
p2s = reduce(vcat,[steadyState[e].nodal.p1[2],
                   steadyState[e].nodal.p2[2]] for e in 1:nElem)
p3s = reduce(vcat,[steadyState[e].nodal.p1[3],
                   steadyState[e].nodal.p2[3]] for e in 1:nElem)                   
cns = [steadyAeroVars[e].aeroCoefficients.cn for e in 1:nElem]

# Create and solve dynamic problem
dynProblem = create_DynamicProblem(model=HealySideslipFFWTOMCGustFloating,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x)
solve!(dynProblem)

# Unpack numerical solution
sol = dynProblem.solutionHistory
states = sol.states
hacStates = sol.hacStates
aeroVars = sol.aeroVariables
t = sol.t
tipOOP = [states[i][nElem].nodal.u2[3] for i in eachindex(t)]
M2_root = [states[i][1].nodal.M1[2] for i in eachindex(t)]
root_αₑ = [rad2deg(aeroVars[i][1].unsteadyKin.αₑ) for i in eachindex(t)]
tqSpan_αₑ = [rad2deg(aeroVars[i][nElementsInner].unsteadyKin.αₑ) for i in eachindex(t)]
tip_αₑ = [rad2deg(aeroVars[i][nElem].unsteadyKin.αₑ) for i in eachindex(t)]
root_cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
tqSpan_cn = [aeroVars[i][nElementsInner].aeroCoefficients.cn for i in eachindex(t)]
tip_cn = [aeroVars[i][nElem].aeroCoefficients.cn for i in eachindex(t)]
root_cm = [aeroVars[i][1].aeroCoefficients.cm for i in eachindex(t)]
tqSpan_cm = [aeroVars[i][nElementsInner].aeroCoefficients.cm for i in eachindex(t)]
tip_cm = [aeroVars[i][nElem].aeroCoefficients.cm for i in eachindex(t)]
root_ct = [aeroVars[i][1].aeroCoefficients.ct for i in eachindex(t)]
tqSpan_ct = [aeroVars[i][nElementsInner].aeroCoefficients.ct for i in eachindex(t)]
tip_ct = [aeroVars[i][nElem].aeroCoefficients.ct for i in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Steady p1
plt_p1 = plot(xlabel="\$x_1/L\$", ylabel="Steady \$p_1\$")
plot!(x1_n/L, p1s, c=:black, lw=lw, label=false)
display(plt_p1)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_steady_p1.pdf"))

# Steady p2
plt_p2 = plot(xlabel="\$x_1/L\$", ylabel="Steady \$p_2\$")
plot!(x1_n/L, p2s, c=:black, lw=lw, label=false)
display(plt_p2)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_steady_p2.pdf"))

# Steady p3
plt_p3 = plot(xlabel="\$x_1/L\$", ylabel="Steady \$p_3\$")
plot!(x1_n/L, p3s, c=:black, lw=lw, label=false)
display(plt_p3)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_steady_p3.pdf"))

# Steady cn over span
plt_cns = plot(xlabel="\$x_1/L\$", ylabel="Steady \$c_n\$")
plot!(x1_e/L, cns, c=:black, lw=lw, label=false)
display(plt_cns)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_steady_cns.pdf"))

# Tip displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [m]")
plot!(t, tipOOP, c=:black, lw=lw, label=false)
display(plt_tipOOP)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_tipOOP.pdf"))

# Root OOP bending moment increment
plt_ΔM2root = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]")
plot!(t, -(M2_root .- M2_root[1]), c=:black, lw=lw, label=false)
display(plt_ΔM2root)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_ΔM2root.pdf"))

# AoA
plt_aoa = plot(xlabel="Time [s]", ylabel="Angle of attack [deg]")
plot!(t, root_αₑ, lw=lw, label="Root")
plot!(t, tqSpan_αₑ, lw=lw, label="3/4-span")
plot!(t, tip_αₑ, lw=lw, label="Tip")
display(plt_aoa)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_aoa.pdf"))

# cn
plt_cn = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, root_cn, lw=lw, label="Root")
plot!(t, tqSpan_cn, lw=lw, label="3/4-span")
plot!(t, tip_cn, lw=lw, label="Tip")
display(plt_cn)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_cn.pdf"))

# cm
plt_cm = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, root_cm, lw=lw, label="Root")
plot!(t, tqSpan_cm, lw=lw, label="3/4-span")
plot!(t, tip_cm, lw=lw, label="Tip")
display(plt_cm)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_cm.pdf"))

# ct
plt_ct = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, root_ct, lw=lw, label="Root")
plot!(t, tqSpan_ct, lw=lw, label="3/4-span")
plot!(t, tip_ct, lw=lw, label="Tip")
display(plt_ct)
savefig(string(absPath,"/HealySideslipFFWTOMCGustFloating_ct.pdf"))

println("Finished HealySideslipFFWTOMCGustFloating.jl")
