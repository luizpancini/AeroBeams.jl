using AeroBeams

# Flag to save figures
saveFigures = true

# Elevator perturbation range
ΔδRange = -π/180*vcat(3:1:6)

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solver
aeroSolver = BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)

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
P = 200

# Airspeed [m/s]
U = 40*0.3048

# Elevator profile variables
tδinit = 1
tδramp = 1
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp

# Time variables
Δt = 5e-2
tf = 120
t_len = round(Int,tf/Δt+1)

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

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(ΔδRange))
dynamicProblem = Array{DynamicProblem}(undef,length(ΔδRange))
t = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
rootPitch = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
rootAoA = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
Δu3 = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
airspeed = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
cn = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
cm = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]
root_κ2 = [fill(NaN, t_len) for _ in 1:length(ΔδRange)]

# Sweep elevator perturbation
for (i,Δδ) in enumerate(ΔδRange)
    println("Solving for Δδ = $(Δδ*180/pi) deg")
    # Model for trim problem
    trimModel,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
    # Create and solve trim problem
    trimProblem[i] = create_TrimProblem(model=trimModel,systemSolver=NRtrim)
    solve!(trimProblem[i])
    # Extract trim variables
    trimAoA = (trimProblem[i].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N)*180/π
    trimThrust = trimProblem[i].x[end-1]*trimProblem[i].model.forceScaling
    trimδ = trimProblem[i].x[end]
    println("Trim variables: AoA = $trimAoA, T = $trimThrust, δ = $(trimδ*180/π)")
    # Set checked elevator deflection profile
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
    dynamicModel,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δ=δ,thrust=trimThrust)
    # Create and solve dynamic problem
    dynamicProblem[i] = create_DynamicProblem(model=dynamicModel,x0=trimProblem[i].x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NRdyn)
    solve!(dynamicProblem[i])
    # Outputs
    Nt = length(dynamicProblem[i].savedTimeVector)
    t[i][1:Nt] = dynamicProblem[i].savedTimeVector
    rootPitch[i][1:Nt] = [dynamicProblem[i].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
    rootAoA[i][1:Nt] = [dynamicProblem[i].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.αₑ-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
    Δu3[i][1:Nt] = [dynamicProblem[i].nodalStatesOverTime[k][midSpanElem].u_n2[3] for k in 1:Nt] .- dynamicProblem[i].nodalStatesOverTime[1][midSpanElem].u_n2[3]
    airspeed[i][1:Nt] = [dynamicProblem[i].aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for k in 1:Nt]
    cn[i][1:Nt] = [dynamicProblem[i].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cn for k in 1:Nt]
    cm[i][1:Nt] = [dynamicProblem[i].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cm for k in 1:Nt]
    root_κ2[i][1:Nt] = [dynamicProblem[i].compElementalStatesOverTime[k][midSpanElem].κ[2] for k in 1:Nt]
end

# Set paths
relPath = "/dev/helios/figures/heliosPitchManeuverDeltaSweepSingle"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
colors = cgrad(:rainbow, length(ΔδRange), categorical=true)
δistr = round(Int,ΔδRange[1]*180/π)
δfstr = round(ΔδRange[end]*180/π)
if typeof(aeroSolver) == Indicial
    aeroSolverStr = "AF"
else
    aeroSolverStr = "DS"
end

# Root angle of attack
plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,tf], ylims=[0,45], yticks=vcat(-90:15:90), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
for (i,Δδ) in enumerate(ΔδRange)
    plot!(t[i], rootAoA[i]*180/π, c=colors[i], lw=lw, label=string("\$\\delta_f\$ = ", round(Int,Δδ*180/π), "\$^\\circ\$"))
end
display(plt_AoA)
savefig(plt_AoA,string(absPath,"/heliosPitchManeuverDeltaSweepSingle_AoAcurves_",aeroSolverStr,"_P",round(P),"_delta",δistr,"_to_",δfstr,".pdf"))

# cm
plt_cma = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs))
for (i,Δδ) in enumerate(ΔδRange)
    plot!(rootPitch[i]*180/π, cm[i], c=colors[i], lw=lw, label=false)
end
display(plt_cma)
savefig(plt_cma,string(absPath,"/heliosPitchManeuverDeltaSweepSingle_cmcurves_",aeroSolverStr,"_P",round(P),"_delta",δistr,"_to_",δfstr,".pdf"))

println("Finished heliosPitchManeuverDeltaSweepSingle.jl")