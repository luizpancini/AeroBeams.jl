using AeroBeams

# Flag to save figures
saveFigures = true

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solvers
aeroSolver = Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction)

# Stiffness factor
λ = 1

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = false

# TF to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Number of elements of beams
nElemStraightSemispan = 10
nElemDihedralSemispan = 5
nElemPod = 2

# Payload [lb]
P = 220

# Airspeed [m/s]
U = 40*0.3048

# Gust variables (duration and intensity as a fraction of forward speed)
τ = 1
γ = 0.01

# DARPA gust
t0 = 1
Ug = U*γ
gustLength = U*τ
L = 0.3048*40*(2+cosd(10))
gustWidth = 2*L
c0 = [0; U*t0; 0]
gust = create_DiscreteSpaceGust(type="DARPA",length=gustLength,width=gustWidth,verticalVelocity=Ug,c0=c0)

# Time variables
Δt = 5e-2
tf = 60

# Set NR system solver for trim problem
relaxFactor = 0.5
maxIter = 100
relTol = 1e-12
absTol = 1e-12
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,relativeTolerance=relTol,absoluteTolerance=absTol)

# Set NR system solver for dynamic problem
allowAdvanceThroughUnconvergedAeroStates = true
maxIter = 50
alwaysUpdateJacobian = false
minConvRateAeroJacUpdate = 1.2
minConvRateJacUpdate = 1.2
rtol = 1e-8
NRdyn = create_NewtonRaphson(maximumIterations=maxIter,alwaysUpdateJacobian=alwaysUpdateJacobian,minConvRateAeroJacUpdate=minConvRateAeroJacUpdate,minConvRateJacUpdate=minConvRateJacUpdate,relativeTolerance=rtol,allowAdvanceThroughUnconvergedAeroStates=allowAdvanceThroughUnconvergedAeroStates,displayStatus=false)

# Model for trim problem
trimModel,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,gust=gust)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=trimModel,systemSolver=NRtrim)
solve!(trimProblem)

# Extract trim variables
trimAoA = (trimProblem.aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N)*180/π
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]
println("Trim variables: AoA = $trimAoA, T = $trimThrust, δ = $(trimδ*180/π)")

# Model for dynamic problem
dynamicModel,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,gust=gust)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=dynamicModel,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NRdyn)
solve!(dynamicProblem)

# Outputs
t = dynamicProblem.savedTimeVector
Nt = length(t)
rootPitch = [dynamicProblem.aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
rootAoA = [dynamicProblem.aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.αₑ-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
airspeed = [dynamicProblem.aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for k in 1:Nt]

# Set paths
relPath = "/dev/helios/figures/heliosDarpaGustSingle"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0

# Root angle of attack
plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
        plot!(t, rootAoA*180/π, c=:black, lw=lw, label=false)
display(plt_AoA)

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [ft/s]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, airspeed, c=:black, lw=lw, label=false)
display(plt_U)

# Strings
τstr = round(τ,digits=1)
γstr = round(Int,γ*100)

# Save figures, if applicable
if saveFigures
    savefig(plt_AoA,string(absPath,"/heliosDarpaGustSingle_AoA_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_U,string(absPath,"/heliosDarpaGustSingle_U_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
end

println("Finished heliosDarpaGustSingle.jl")