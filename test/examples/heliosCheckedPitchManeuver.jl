using AeroBeams

# Stiffness factor
λ = 1

# Wing airfoil
wingAirfoil = NACA23012A

# Option for reduced chord
reducedChord = false

# TF to include beam pods and number of elements
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Number of elements of beams
nElemStraightSemispan = 10
nElemPod = 2

# Aerodynamic solver
aeroSolver = Indicial()

# Payload [lb]
P = 200

# Airspeed [m/s]
U = 40*0.3048

# Set NR system solver for trim problem
NR = create_NewtonRaphson(ρ=0.5,pseudoInverseMethod=:MoorePenrose,relativeTolerance=1e-12,maximumIterations=100,displayStatus=false)

# Model for trim problem
heliosTrim,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,nElemStraightSemispan=nElemStraightSemispan,nElemPod=nElemPod,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR)
solve!(trimProblem)

# Extract trim variables
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]
println("Trim variables: T = $(trimThrust), δ = $(trimδ*180/π)")

# Set checked elevator deflection profile
Δδ = -10*π/180
tδinit = 1
tδramp = 1
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp
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
heliosDynamic,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,nElemStraightSemispan=nElemStraightSemispan,nElemPod=nElemPod,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=δ,thrust=trimThrust,reducedChord=reducedChord,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing)

# Time variables
Δt = 5e-2
tf = 10

# Set NR system solver for dynamic problem
maxit = 50
NR = create_NewtonRaphson(maximumIterations=maxit,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=heliosDynamic,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NR)
solve!(dynamicProblem)

# Unpack numerical solution
t = dynamicProblem.timeVector
rootAoA = [dynamicProblem.aeroVariablesOverTime[i][midSpanElem].flowAnglesAndRates.αₑ for i in 1:length(t)]
Δu3 = [dynamicProblem.nodalStatesOverTime[i][midSpanElem].u_n2[3] for i in 1:length(t)] .- dynamicProblem.nodalStatesOverTime[1][midSpanElem].u_n2[3]

println("Finished heliosCheckedPitchManeuver.jl")