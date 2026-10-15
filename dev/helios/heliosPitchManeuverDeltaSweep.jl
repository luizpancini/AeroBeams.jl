using AeroBeams

# Elevator perturbation range
ΔδRange = π/180*vcat(5:5:30)

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solvers
aeroSolvers = [Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction); BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)]

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
P = 0

# Airspeed [m/s]
U = 40*0.3048

# Elevator profile variables
tδinit = 1
tδramp = 1
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp

# Time variables
Δt = 5e-2
tf = 60
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
trimProblem = Array{TrimProblem}(undef,length(aeroSolvers),length(ΔδRange))
dynamicProblem = Array{DynamicProblem}(undef,length(aeroSolvers),length(ΔδRange))
t = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
rootPitch = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
rootAoA = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
Δu3 = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
airspeed = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
cn = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
cm = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]
root_κ2 = [fill(NaN, t_len) for _ in 1:length(aeroSolvers), _ in 1:length(ΔδRange)]

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep elevator perturbation
    for (j,Δδ) in enumerate(ΔδRange)
        println("Solving for aeroSolver $i, Δδ = $(Δδ*180/pi) deg")
        # Model for trim problem
        trimModel,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=trimModel,systemSolver=NRtrim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA = (trimProblem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N)*180/π
        trimThrust = trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling
        trimδ = trimProblem[i,j].x[end]
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
        dynamicProblem[i,j] = create_DynamicProblem(model=dynamicModel,x0=trimProblem[i,j].x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NRdyn)
        solve!(dynamicProblem[i,j])
        # Outputs
        Nt = length(dynamicProblem[i,j].savedTimeVector)
        t[i,j][1:Nt] = dynamicProblem[i,j].savedTimeVector
        rootPitch[i,j][1:Nt] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
        rootAoA[i,j][1:Nt] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.αₑ-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
        Δu3[i,j][1:Nt] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[3]
        airspeed[i,j][1:Nt] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for k in 1:Nt]
        cn[i,j][1:Nt] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cn for k in 1:Nt]
        cm[i,j][1:Nt] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cm for k in 1:Nt]
        root_κ2[i,j][1:Nt] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[2] for k in 1:Nt]
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosPitchManeuverDeltaSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0

δistr = round(Int,ΔδRange[1])
δfstr = round(Int,ΔδRange[end])

# Mesh grid
T = t[1]' .* ones(length(ΔδRange))
D = ones(length(t[1]))' .* ΔδRange*180/pi
Z = 180/pi*hcat(rootAoA[1,:]...)'

# Root angle of attack surface
plotlyjs()
surf_AoA = surface(xlabel="Time [s]", ylabel="δ [deg]", zlabel="Root AoA [deg]", xlims=[0,tf], ylims=[0,30], zlims=[-10,90], tickfont=font(ts), guidefont=font(fs), colorbar=false, yticks=vcat(0:10:30), zticks=vcat(0:10:90))
surface!(T, D, Z)
display(surf_AoA)
savefig(surf_AoA,string(absPath,"/heliosPitchManeuverDeltaSweep_AoAsurf_P",round(P),"_delta",δistr,"_to_",δfstr,".pdf"))

# Root angle of attack of selected perturbations
gr()
Δδ2plot = [5,10,15,20,25,30]
colors = cgrad(:rainbow, length(Δδ2plot), categorical=true)
solverLabels = ["_attached" "_ds"]
for (i,aeroSolver) in enumerate(aeroSolvers)
    n = 0
    plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,tf], ylims=[-15,90], yticks=vcat(-15:15:90), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
    for (j,Δδ) in enumerate(ΔδRange)
        if !(round(Int,Δδ*180/π) in Δδ2plot)
            continue
        else
            n+=1
        end
        if i==1
            plot!(t[i,j], rootAoA[i,j]*180/π, c=colors[n], lw=lw, label=string("\$\\delta_f\$ = ", round(Int,Δδ*180/π), "\$^\\circ\$"))
        else
            plot!(t[i,j], rootAoA[i,j]*180/π, c=colors[n], lw=lw, label=false)
        end
    end
    display(plt_AoA)
    savefig(plt_AoA,string(absPath,"/heliosPitchManeuverDeltaSweep_AoAcurves_P",round(P),"_delta",δistr,"_to_",δfstr,"_",solverLabels[i],".pdf"))
end

println("Finished heliosPitchManeuverDeltaSweep.jl")