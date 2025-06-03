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
P = 150

# Airspeed [m/s]
U = 40*0.3048

# Gust variables (duration and intensity as a fraction of forward speed)
τ = 27.5
γ = 0.1

# DARPA gust
t0 = 1
Ug = U*γ
gustLength = U*τ
L = 0.3048*40*(2+cosd(10))
gustWidth = 2*L
c0 = [0; U*t0; 0]
gust = create_DiscreteSpaceGust(type="DARPA",length=gustLength,width=gustWidth,verticalVelocity=Ug,c0=c0)

# Time variables
Δt = [5e-2 5e-2; 5e-2 5e-2]
tf = [90 90; 90 90]

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

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep stiffness factors
    for (j,λ) in enumerate(λRange)
        println("Solving for aeroSolver $i, λ = $λ")
        # Model for trim problem
        trimModel,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,gust=gust)
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=trimModel,systemSolver=NRtrim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA = (trimProblem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N)*180/π
        trimThrust = trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling
        trimδ = trimProblem[i,j].x[end]
        println("Trim variables: AoA = $trimAoA, T = $trimThrust, δ = $(trimδ*180/π)")
        # Model for dynamic problem
        dynamicModel,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,gust=gust)
        # Create and solve dynamic problem
        dynamicProblem[i,j] = create_DynamicProblem(model=dynamicModel,x0=trimProblem[i,j].x[1:end-2],finalTime=tf[i,j],Δt=Δt[i,j],skipInitialStatesUpdate=true,systemSolver=NRdyn)
        solve!(dynamicProblem[i,j])
        # Outputs
        t[i,j] = dynamicProblem[i,j].savedTimeVector
        Nt = length(t[i,j])
        rootPitch[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
        rootAoA[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.αₑ-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
        Δu3[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[3]
        airspeed[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for k in 1:Nt]
        cn[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cn for k in 1:Nt]
        cm[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cm for k in 1:Nt]
        root_κ2[i,j] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[2] for k in 1:Nt]
        tip_Δu3[i,j] = [(dynamicProblem[i,j].nodalStatesOverTime[k][1].u_n1[3] - dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3]) for k in 1:Nt]
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosDarpaGust"
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
labels = ["Attached flow - Elastic" "Attached flow - Rigid"; "Dynamic stall - Elastic" "Dynamic stall - Rigid"]
ls = [:solid, :dash]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)

# Root angle of attack
plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,maximum(tf)], ylims=[-15,90], yticks=vcat(-15:15:90), tickfont=font(ts), guidefont=font(fs), legend=(0.6,0.9), legendfontsize=lfs)
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], rootAoA[i,j]*180/π, c=colors[i], ls=ls[j], lw=lw, label=labels[i,j])
    end
end
display(plt_AoA)

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [ft/s]", xlims=[0,maximum(tf)], ylims=[0,Inf], yticks=vcat(0:20:1e3), tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], airspeed[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_U)

# Elevation change
plt_deltaz = plot(xlabel="Time [s]", ylabel="Elevation change [m]", xlims=[0,maximum(tf)], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], Δu3[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_deltaz)

# Root bending curvature
plt_k2 = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", xlims=[0,maximum(tf)], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], root_κ2[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_k2)

# Tip OOP displacement
L = 0.3048*40*(2+cosd(10))
plt_OOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,maximum(tf)], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], tip_Δu3[i,j]/L*100, c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_OOP)

# cn x time
plt_cnt = plot(xlabel="Time [s]", ylabel="\$c_n\$", xlims=[0,maximum(tf)], ylims=[-0.5,3], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], cn[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cnt)

# cm  x time
plt_cmt = plot(xlabel="Time [s]", ylabel="\$c_m\$", xlims=[0,maximum(tf)], ylims=[-0.6,0.2], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(t[i,j], cm[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cmt)

# cn x α
plt_cna = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_n\$", xlims=[-Inf,Inf], ylims=[-Inf,Inf], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(rootPitch[i,j]*180/π, cn[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cna)

# cm  x α
plt_cma = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_m\$", xlims=[-Inf,Inf], ylims=[-Inf,Inf], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(rootPitch[i,j]*180/π, cm[i,j], c=colors[i], ls=ls[j], lw=lw, label=false)
    end
end
display(plt_cma)

# Strings
τstr = round(τ,digits=1)
γstr = round(Int,γ*100)

# Save figures, if applicable
if saveFigures
    savefig(plt_AoA,string(absPath,"/heliosDarpaGust_AoA_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_U,string(absPath,"/heliosDarpaGust_U_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_deltaz,string(absPath,"/heliosDarpaGust_deltaz_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_k2,string(absPath,"/heliosDarpaGust_k2_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_OOP,string(absPath,"/heliosDarpaGust_OOP_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_cnt,string(absPath,"/heliosDarpaGust_cnt_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_cmt,string(absPath,"/heliosDarpaGust_cmt_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_cna,string(absPath,"/heliosDarpaGust_cna_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    savefig(plt_cma,string(absPath,"/heliosDarpaGust_cma_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
end

# Animations
plot_dynamic_deformation(dynamicProblem[1,1],refBasis="I",view=(60,15),plotDistLoads=false,plotFrequency=2,plotLimits=([-40,40],[-5,1000],[-40,10]),save=saveFigures,savePath=string(relPath,"/heliosDarpaGust_attached_elastic_P",round(P),"_tau",τstr,"_gamma",γstr,".gif"),displayProgress=true)
plot_dynamic_deformation(dynamicProblem[2,1],refBasis="I",view=(60,15),plotDistLoads=false,plotFrequency=2,plotLimits=([-40,40],[-5,800],[-200,10]),save=saveFigures,savePath=string(relPath,"/heliosDarpaGust_ds_elastic_P",round(P),"_tau",τstr,"_gamma",γstr,".gif"),displayProgress=true)

println("Finished heliosDarpaGust.jl")