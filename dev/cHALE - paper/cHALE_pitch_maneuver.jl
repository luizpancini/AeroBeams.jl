using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Altitude [m]
h = 20e3

# Airspeed [m/s]
U = 25

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Stiffness factor (for the structure)
λ = 1

# Pre-curvatures
k1 = 0
k2 = 0.045

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
σstep = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)

# Model for trim problem
cHALEtrim,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k1=k1,k2=k2,hasInducedDrag=hasInducedDrag)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR)
solve!(trimProblem)

# Extract trim variables and outputs
trimAoA = (trimProblem.aeroVariablesOverσ[end][cHALEtrim.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblem.aeroVariablesOverσ[end][cHALEtrim.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]
println("Trim outputs: AoA = $(trimAoA*180/π), T = $(trimThrust), δ = $(trimδ*180/π)")

# Set checked elevator deflection profile
Δδ = -10*π/180
tδinit = 0.5
tδramp = 0.5
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
cHALEdynamic,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,k1=k1,k2=k2,δElev=δ,thrust=trimThrust)

# Time variables
Δt = 1e-2
tf = 120

# Set NR system solver for dynamic problem
maxit = 100
NR = create_NewtonRaphson(maximumIterations=maxit)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=cHALEdynamic,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NR)
solve!(dynamicProblem)

# Get wing root elements
lRootElem = div(nElemWing,2)
rRootElem = lRootElem+1

# Unpack numerical solution
t = dynamicProblem.timeVector
wingAoA = [(dynamicProblem.aeroVariablesOverTime[i][lRootElem].flowAnglesAndRates.αₑ+dynamicProblem.aeroVariablesOverTime[i][rRootElem].flowAnglesAndRates.αₑ)/2 for i in 1:length(t)]
airspeed = [(dynamicProblem.aeroVariablesOverTime[i][lRootElem].flowVelocitiesAndRates.U∞+dynamicProblem.aeroVariablesOverTime[i][rRootElem].flowVelocitiesAndRates.U∞)/2 for i in 1:length(t)]

# Set paths
relPathFig = "/dev/outputs/figures/cHALE_pitch_maneuver"
relPathData = "/dev/outputs/data/cHALE_pitch_maneuver"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
mkpath(absPathFig)
mkpath(absPathData)

# Plots
using Plots, DelimitedFiles
lw = 2
gr()

# Body AoA
plt1 = plot(xlabel="Time [s]", ylabel="Normalized wing root angle of attack")
plot!(t, wingAoA./wingAoA[1], c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_AoA_U",U,"_k2",k2,".pdf")))

# Airspeed
plt2 = plot(xlabel="Time [s]", ylabel="Normalized airspeed")
plot!(t, airspeed./airspeed[1], c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_airspeed_U",U,"_k2",k2,".pdf")))

# Save arrays
writedlm(string(absPathData,"/cHALE_pitch_maneuver_t_U",U,"_k2",k2,".txt"), t)
writedlm(string(absPathData,"/cHALE_pitch_maneuver_AoA_U",U,"_k2",k2,".txt"), wingAoA)
writedlm(string(absPathData,"/cHALE_pitch_maneuver_airspeed_U",U,"_k2",k2,".txt"), airspeed)

println("Finished cHALE_pitch_maneuver.jl")