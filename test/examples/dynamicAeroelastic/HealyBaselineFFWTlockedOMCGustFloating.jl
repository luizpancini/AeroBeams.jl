using AeroBeams

# Gust frequency range [Hz]
ωRange = [1,3,10]

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Hinge configuration
hingeConfiguration = :locked

# Fold angle [rad]
foldAngle = nothing

# Root pitch angle [rad]
θ = deg2rad(7.5)

# Airspeed [m/s]
U = 18

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 12

# Time variables
Δt = 2.5e-3
tf = 2

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Initialize outputs
t = Array{Vector{Float64}}(undef,length(ωRange))
M2root = Array{Vector{Float64}}(undef,length(ωRange))

# Loop gust frequency range
for (i,ω) in enumerate(ωRange)
    # Display progress
    display("Solving for ω = $ω Hz")
    # Gust
    gust = create_OneMinusCosineGust(duration=1/ω,verticalVelocity=Ug)
    # Model
    model = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,foldAngle=foldAngle,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust) 
    # Solve steady problem for initial conditions
    steadyProblem = create_SteadyProblem(model=model,systemSolver=systemSolver)
    solve!(steadyProblem)
    # Create and solve dynamic problem
    dynProblem = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x)
    solve!(dynProblem)
    # Unpack numerical solution
    sol = dynProblem.solutionHistory
    states = sol.states
    t[i] = sol.t
    M2root[i] = [states[j][1].nodal.M1[2] for j in eachindex(t[i])]
end

# Load reference data
ΔWRBM_Healy_1Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedOMCGustFloating/DWRBM_baseline_1Hz.txt")
ΔWRBM_Healy_3Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedOMCGustFloating/DWRBM_baseline_3Hz.txt")
ΔWRBM_Healy_10Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedOMCGustFloating/DWRBM_baseline_10Hz.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTlockedOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(ωRange), categorical=true)
ts = 10
fs = 16
lw = 2
gr()

# Root OOP bending moment increment
plt_ΔWRBM = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
for (i,ω) in enumerate(ωRange)
    plot!(t[i], -(M2root[i] .- M2root[i][1]), c=colors[i], lw=lw, ls=:solid, label="\$\\omega_g = $ω\$ Hz")
    if ω==1
        plot!(ΔWRBM_Healy_1Hz[1,:],ΔWRBM_Healy_1Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif ω==3
        plot!(ΔWRBM_Healy_3Hz[1,:],ΔWRBM_Healy_3Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif ω==10
        plot!(ΔWRBM_Healy_10Hz[1,:],ΔWRBM_Healy_10Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    end
end
display(plt_ΔWRBM)
savefig(string(absPath,"/HealyBaselineFFWTlockedOMCGustFloating_DWRBM.pdf"))

println("Finished HealyBaselineFFWTlockedOMCGustFloating.jl")