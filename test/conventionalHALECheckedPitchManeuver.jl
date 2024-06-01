using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Airspeed [m/s]
U = 25

# Wing and stabilizers parasite drag
wingCd0 = stabsCd0 = 1e-2

# Discretization
nElemWing = 20

# Set NR system solver for trim problem
NR = create_NewtonRaphson(ρ=0.5,relativeTolerance=1e-12,maximumIterations=50,displayStatus=false)

# Model for trim problem
conventionalHALEtrim,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElemWing=nElemWing,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=conventionalHALEtrim,systemSolver=NR)
solve!(trimProblem)

# Extract trim variables
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]

# Set checked elevator deflection profile
Δδ = 5*π/180
tδinit = 0.5
tδpeak = 1+tδinit
tδfinal = 1+tδpeak
δ = t -> ifelse(
    t <= tδinit, 
    trimδ,
    ifelse(
        t <= tδpeak, 
        trimδ + Δδ * ((t-tδinit) / (tδpeak-tδinit)),
        ifelse(
            t <= tδfinal, 
            trimδ + Δδ - Δδ * ((t-tδpeak) / (tδfinal-tδpeak)),
            trimδ
        )
    )
)

# Model for dynamic problem
conventionalHALEdynamic,leftWing,rightWing,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElemWing=nElemWing,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=δ,thrust=trimThrust)

# Time variables
Δt = 1e-3
tf = 5

# Set NR system solver for trim problem
maxit = 100
NR = create_NewtonRaphson(maximumIterations=maxit,displayStatus=false)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=conventionalHALEdynamic,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NR)
solve!(dynamicProblem)
# @profview solve!(dynamicProblem)

# Get wing root elements
lRootElem = div(nElemWing,2)
rRootElem = lRootElem+1

# Unpack numerical solution
t = dynamicProblem.timeVector
rootAoA = [(dynamicProblem.flowVariablesOverTime[i][lRootElem].αₑ+dynamicProblem.flowVariablesOverTime[i][rRootElem].αₑ)/2 for i in 1:length(t)]
Δu3 = [dynamicProblem.nodalStatesOverTime[i][lRootElem].u_n2[3] for i in 1:length(t)] .- dynamicProblem.nodalStatesOverTime[1][lRootElem].u_n2[3]

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALECheckedPitchManeuver_altitude.pdf"))
# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, rootAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALECheckedPitchManeuver_rootAoA.pdf"))

println("Finished conventionalHALECheckedPitchManeuver.jl")