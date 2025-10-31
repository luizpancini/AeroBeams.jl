using AeroBeams

# Gust frequency [Hz]
ω = 3

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Hinge configuration
hingeConfiguration = "free"

# Root pitch angle [rad]
θ = 7.5*π/180

# Airspeed [m/s]
U = 18

# Stiffness of the spring around the hinge
kSpring = 1e-4
kIPBendingHinge = 1e-1

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 12

# Solution method for hinge constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = true

# Time variables
Δt = 2e-3
tf = 1

# System solver
σ0 = 1
maxIter = 200
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Gust
gust = create_OneMinusCosineGust(initialTime=0,duration=1/ω,verticalVelocity=Ug)

# Model
model = create_HealyBaselineFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,hingeConfiguration=hingeConfiguration,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust)

# Solve steady problem for initial conditions
steadyProblem = create_SteadyProblem(model=model,systemSolver=NR)
solve!(steadyProblem)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x)
solve!(problem)

# Unpack numerical solution
t = problem.savedTimeVector
M2root = [problem.nodalStatesOverTime[j][1].M_n1[2] for j in 1:length(t)]
ϕ = [problem.hingeAxisConstraintsDataOverTime[j][1].ϕ*180/π for j in 1:length(t)]

# Do plots if running locally
is_ci = get(ENV, "CI", "false") == "true" || get(ENV, "GITHUB_ACTIONS", "false") == "true"
if !is_ci

    using Plots, ColorSchemes

    # Plot configurations
    ts = 10
    fs = 16
    lfs = 12
    lw = 2
    gr()

    # Fold angle increment
    plt_ϕ = plot(xlabel="Time [s]", ylabel="Fold angle increment [deg]", ylims=[-20,20], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
    plot!(t, -(ϕ .- ϕ[1]), c=:black, lw=lw, ls=:solid, label=false)
    display(plt_ϕ)

    # Root OOP bending moment increment
    plt_ΔWRBM = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", tickfont=font(ts), guidefont=font(fs))
    plot!(t, -(M2root .- M2root[1]), c=:black, lw=lw, ls=:solid, label=false)
    display(plt_ΔWRBM)
end

println("Finished HealyBaselineFFWTOMCGustFloating2.jl")