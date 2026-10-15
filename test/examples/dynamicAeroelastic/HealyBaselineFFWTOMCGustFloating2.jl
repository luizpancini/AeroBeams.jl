using AeroBeams

# Gust frequency [Hz]
ω = 3

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Hinge configuration
hingeConfiguration = :free

# Root pitch angle [rad]
θ = deg2rad(7.5)

# Airspeed [m/s]
U = 18

# Discretization
nElementsInner = 16
nElementsFFWT = 8
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 17

# Time variables
Δt = 2.5e-3
tf = 1

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Gust
gust = create_OneMinusCosineGust(duration=1/ω,verticalVelocity=Ug)

# Model
model = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust)

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=model,systemSolver=systemSolver)
solve!(steadyProblem)

# Create and solve dynamic problem
dynProblem = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x)
solve!(dynProblem)

# Unpack numerical solution
sol = dynProblem.solutionHistory
states = sol.states
hacStates = sol.hacStates
t = sol.t
M2root = [states[j][1].nodal.M1[2] for j in eachindex(t)]
ϕ = [rad2deg(hacStates[j][1].ϕ) for j in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTOMCGustFloating2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lfs = 12
lw = 2

# Fold angle increment
plt_ϕ = plot(xlabel="Time [s]", ylabel="Fold angle increment [deg]", ylims=[-20,20], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
plot!(t, -(ϕ .- ϕ[1]), c=:black, lw=lw, ls=:solid, label=false)
display(plt_ϕ)
savefig(string(absPath,"/HealyBaselineFFWTOMCGustFloating2_fold_angle_increment.pdf"))

# Root OOP bending moment increment
plt_ΔWRBM = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", tickfont=font(ts), guidefont=font(fs))
plot!(t, -(M2root .- M2root[1]), c=:black, lw=lw, ls=:solid, label=false)
display(plt_ΔWRBM)
savefig(string(absPath,"/HealyBaselineFFWTOMCGustFloating2_DWRBM.pdf"))

println("Finished HealyBaselineFFWTOMCGustFloating2.jl")