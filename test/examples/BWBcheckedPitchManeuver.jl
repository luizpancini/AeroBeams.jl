using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Flight altitude [m]
h = 0e3

# Airspeed [m/s]
U = 80

# Set NR system solver for trim problem
NRtrim = create_NewtonRaphson(ρ=0.5,relativeTolerance=1e-12,maximumIterations=100,displayStatus=false)

# Model for trim problem
BWBtrim = create_BWB(aeroSolver=aeroSolver,δElevIsTrimVariable=true,thrustIsTrimVariable=true,altitude=h,airspeed=U)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=BWBtrim,systemSolver=NRtrim)
solve!(trimProblem)

# Extract trim variables
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]
println("Trim variables: T = $(trimThrust), δ = $(trimδ*180/π)")

# Set checked elevator deflection profile
Δδ = -2*π/180
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
BWBdynamic = create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElev=δ,thrust=trimThrust)

# Time variables
Δt = 5e-3
tf = 5

# Set NR system solver for dynamic problem
maxit = 100
NRdyn = create_NewtonRaphson(maximumIterations=maxit,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=BWBdynamic,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NRdyn)
solve!(dynamicProblem)

# Unpack numerical solution
t = dynamicProblem.timeVector
rootAoA = [(dynamicProblem.aeroVariablesOverTime[i][11].flowAnglesAndRates.αₑ+dynamicProblem.aeroVariablesOverTime[i][12].flowAnglesAndRates.αₑ)/2 for i in 1:length(t)]
Δu3 = [dynamicProblem.nodalStatesOverTime[i][11].u_n2[3] for i in 1:length(t)] .- dynamicProblem.nodalStatesOverTime[1][11].u_n2[3]

println("Finished BWBcheckedPitchManeuver.jl")