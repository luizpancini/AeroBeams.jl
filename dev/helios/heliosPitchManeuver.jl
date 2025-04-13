using AeroBeams

# Flag to save figures
saveFigures = true

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solvers
aeroSolvers = [Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction); BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)]

# Stiffness factor range
λRange = [1, 1e4]

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = true

# TF to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Number of elements of beams
nElemStraightSemispan = 10
nElemDihedralSemispan = 5
nElemPod = 2

# Payload [lb]
P = 100

# Airspeed [m/s]
U = 40*0.3048

# Elevator profile variables
Δδ = 8*π/180
tδinit = 1
tδramp = 1
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp

# Time variables
Δt = [5e-2 1e-2; 5e-2 1e-2]
tf = [120 120; 30 120]

# Set NR system solver for trim problem
relaxFactor = 0.5
maxIter = 50
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=false)

# Set NR system solver for dynamic problem
allowAdvanceThroughUnconvergedAeroStates = true
maxIter = 50
alwaysUpdateJacobian = false
minConvRateAeroJacUpdate = 1.2
minConvRateJacUpdate = 1.2
rtol = 1e-8
NRdyn = create_NewtonRaphson(maximumIterations=maxIter,alwaysUpdateJacobian=alwaysUpdateJacobian,minConvRateAeroJacUpdate=minConvRateAeroJacUpdate,minConvRateJacUpdate=minConvRateJacUpdate,relativeTolerance=rtol,allowAdvanceThroughUnconvergedAeroStates=allowAdvanceThroughUnconvergedAeroStates,displayStatus=false)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(aeroSolvers),length(λRange))
dynamicProblem = Array{DynamicProblem}(undef,length(aeroSolvers),length(λRange))
t = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
rootPitch = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
rootAoA = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
Δu3 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
airspeed = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
cn = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
cm = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
root_κ2 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
tip_Δu3 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange))
χ = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers),length(λRange))
χdot = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers),length(λRange))

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep stiffness factors
    for (j,λ) in enumerate(λRange)
        println("Solving for aeroSolver $i, λ = $λ")
        # Model for trim problem
        trimModel,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=trimModel,systemSolver=NRtrim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA = trimProblem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
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
        dynamicProblem[i,j] = create_DynamicProblem(model=dynamicModel,x0=trimProblem[i,j].x[1:end-2],finalTime=tf[i,j],Δt=Δt[i,j],skipInitialStatesUpdate=true,systemSolver=NRdyn)
        solve!(dynamicProblem[i,j])
        # Outputs
        t[i,j] = dynamicProblem[i,j].savedTimeVector
        Nt = length(t[i,j])
        rootPitch[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.α for k in 1:Nt]
        rootAoA[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.αₑ for k in 1:Nt]
        Δu3[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[3]
        airspeed[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞ for k in 1:Nt]
        cn[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cn for k in 1:Nt]
        cm[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cm for k in 1:Nt]
        root_κ2[i,j] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[2] for k in 1:Nt]
        tip_Δu3[i,j] = [(dynamicProblem[i,j].nodalStatesOverTime[k][1].u_n1[3] - dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3]) for k in 1:Nt]
        χ[i,j] = [dynamicProblem[i,j].elementalStatesOverTime[k][midSpanElem].χ for k in 1:Nt]
        χdot[i,j] = [dynamicProblem[i,j].elementalStatesRatesOverTime[k][midSpanElem].χdot for k in 1:Nt]
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosPitchManeuver"
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
labels = ["Linear aero - Elastic" "Linear aero - Rigid"; "Dynamic stall - Elastic" "Dynamic stall - Rigid"]
ls = [:solid, :dash]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)
stateColors = cgrad(:rainbow, 8, categorical=true)

# Root angle of attack
plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,90], ylims=[0,15], tickfont=font(ts), guidefont=font(fs), legend=:topright, legendfontsize=lfs, legend_column=1)
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], rootAoA[i,j]*180/π, c=colors[i], ls=ls[j], lw=lw, label=labels[i,j])
    end
end
display(plt_AoA)

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", xlims=[0,90], ylims=[0,25], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], airspeed[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_U)

# Elevation change
plt_deltaz = plot(xlabel="Time [s]", ylabel="Elevation change [m]", xlims=[0,90], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], Δu3[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_deltaz)

# Root bending curvature
plt_k2 = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", xlims=[0,90], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], root_κ2[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_k2)

# Tip vertical displacement
L = 0.3048*40*(2+cosd(10))
plt_tipdisp = plot(xlabel="Time [s]", ylabel="Tip vert. disp. [% semispan]", xlims=[0,90], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], tip_Δu3[i,j]/L*100, c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_tipdisp)

# cn x time
plt_cnt = plot(xlabel="Time [s]", ylabel="\$c_n\$", xlims=[0,90], ylims=[-0.5,3], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], cn[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cnt)

# cm  x time
plt_cmt = plot(xlabel="Time [s]", ylabel="\$c_m\$", xlims=[0,90], ylims=[-0.6,0.2], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], cm[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cmt)

# cn x α
plt_cna = plot(xlabel="Root pitch angle [deg]", ylabel="\$c_n\$", xlims=[-15,40], ylims=[-0.5,3], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(rootPitch[i,j]*180/π, cn[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cna)

# cm  x α
plt_cma = plot(xlabel="Root pitch angle [deg]", ylabel="\$c_m\$", xlims=[-15,40], ylims=[-0.6,0.3], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(rootPitch[i,j]*180/π, cm[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cma)

# # States for dynamic stall elastic case
# plt_χ = plot(xlabel="Time [s]", ylabel="States", tickfont=font(ts), guidefont=font(fs))
# for (i,aeroSolver) in enumerate(aeroSolvers)
#     for m in eachindex(χ[i,1][1])
#         χ_now = [χ[i,1][k][m] for k in eachindex(t[i,1])]
#         plot!(t[i,1], χ_now, c=stateColors[m], lw=lw, label=false)
#     end
# end
# display(plt_χ)

# # States' rates for dynamic stall elastic case
# plt_χdot = plot(xlabel="Time [s]", ylabel="States rates", tickfont=font(ts), guidefont=font(fs))
# for (i,aeroSolver) in enumerate(aeroSolvers)
#     for m in eachindex(χdot[i,1][1])
#         χdot_now = [χdot[i,1][k][m] for k in eachindex(t[i,1])]
#         plot!(t[i,1], χdot_now, c=stateColors[m], lw=lw, label=false)
#     end
# end
# display(plt_χdot)

# Snapshots
plt_snap = plot_snapshots(dynamicProblem[1,1],refBasis="I",scale=1,plotBCs=false,plotDistLoads=false,plotAxes=true,plotGrid=true,view=(120,15),snapshots=[0,5,10,20,30,60],plotLimits=([-50,50],[-5,800],[-50,5]),save=true,savePath=string(relPath,"/heliosPitchManeuver_snapshots.pdf"))
display(plt_snap)
gr()

# Save figures, if applicable
if saveFigures
    savefig(plt_AoA,string(absPath,"/heliosPitchManeuver_AoA_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_U,string(absPath,"/heliosPitchManeuver_U_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_deltaz,string(absPath,"/heliosPitchManeuver_deltaz_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_k2,string(absPath,"/heliosPitchManeuver_k2_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_tipdisp,string(absPath,"/heliosPitchManeuver_OOP_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_cnt,string(absPath,"/heliosPitchManeuver_cnt_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_cmt,string(absPath,"/heliosPitchManeuver_cmt_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_cna,string(absPath,"/heliosPitchManeuver_cna_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
    savefig(plt_cma,string(absPath,"/heliosPitchManeuver_cma_P",round(P),"_delta",round(Int,Δδ*180/pi),".pdf"))
end

println("Finished heliosPitchManeuver.jl")