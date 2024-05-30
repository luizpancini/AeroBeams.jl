using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Option for reduced chord
reducedChord = true

# Stiffness factor
λ = 1

# TF to include beam pods and number of elements
beamPods = false

# Number of elements of beams
nElemStraightSemispan = 10
nElemPod = 3

# Aerodynamic solver
aeroSolver = Indicial()

# Payload [lb]
P = 0

# Airspeed [m/s]
U = 40*0.3048

# Set NR system solver for trim problem
NR = create_NewtonRaphson(ρ=0.5,relativeTolerance=1e-12,maximumIterations=50,displayStatus=false)

# Model for trim problem
heliosTrim,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,nElemStraightSemispan=nElemStraightSemispan,nElemPod=nElemPod,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR)
solve!(trimProblem)

# Extract trim variables
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]

# Set checked elevator deflection profile
Δδ = 1*π/180
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
heliosDynamic,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,nElemStraightSemispan=nElemStraightSemispan,nElemPod=nElemPod,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=δ,thrust=trimThrust,reducedChord=reducedChord)

# Time variables
Δt = 1e-3
tf = 5

# Set NR system solver for trim problem
maxit = 100
NR = create_NewtonRaphson(maximumIterations=maxit,displayStatus=false)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=heliosDynamic,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NR)
solve!(dynamicProblem)
# @profview solve!(dynamicProblem)

# Unpack numerical solution
t = dynamicProblem.timeVector
rootAoA = [dynamicProblem.flowVariablesOverTime[i][midSpanElem].αₑ for i in 1:length(t)]
Δu3 = [dynamicProblem.nodalStatesOverTime[i][midSpanElem].u_n2[3] for i in 1:length(t)] .- dynamicProblem.nodalStatesOverTime[1][midSpanElem].u_n2[3]

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/heliosCheckedPitchManeuver_altitude.pdf"))
# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, rootAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/heliosCheckedPitchManeuver_rootAoA.pdf"))

println("Finished heliosCheckedPitchManeuver.jl")