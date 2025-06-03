using AeroBeams

# Flags to save figures and create animations
saveFigures = true
createAnimations = false

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solvers
aeroSolvers = [Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction); BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)]

# Stiffness factor range
λRange = [1, 1e4]

# Gust duration range
τRange = vcat(2.5:0.5:5)

# Gust intensity (as a fraction of forward speed)
γ = 0.25

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

# Fixed DARPA gust variables
t0 = 1
Ug = U*γ
L = 0.3048*40*(2+cosd(10))
gustWidth = 2*L
c0 = [0; U*t0; 0]

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
trimProblem = Array{TrimProblem}(undef,length(aeroSolvers),length(λRange),length(τRange))
dynamicProblem = Array{DynamicProblem}(undef,length(aeroSolvers),length(λRange),length(τRange))
t = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
rootPitch = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
rootAoA = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
Δu2 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
Δu3 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
airspeed = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
cn = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
cm = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
root_κ1 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
root_κ2 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))
tip_Δu3 = Array{Vector{Float64}}(undef,length(aeroSolvers),length(λRange),length(τRange))

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep stiffness factors
    for (j,λ) in enumerate(λRange)
        # Sweep gust duration
        for (k,τ) in enumerate(τRange)
            println("Solving for aeroSolver $i, λ = $λ, τ = $τ s")
            # DARPA gust
            gust = create_DiscreteSpaceGust(type="DARPA",length=U*τ,width=gustWidth,verticalVelocity=Ug,c0=c0)
            # Model for trim problem
            trimModel,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,gust=gust)
            # Create and solve trim problem
            trimProblem[i,j,k] = create_TrimProblem(model=trimModel,systemSolver=NRtrim)
            solve!(trimProblem[i,j,k])
            # Extract trim variables
            trimAoA = (trimProblem[i,j,k].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N)*180/π
            trimThrust = trimProblem[i,j,k].x[end-1]*trimProblem[i,j,k].model.forceScaling
            trimδ = trimProblem[i,j,k].x[end]
            println("Trim variables: AoA = $trimAoA, T = $trimThrust, δ = $(trimδ*180/π)")
            # Model for dynamic problem
            dynamicModel,_ = create_Helios(aeroSolver=aeroSolver,reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,gust=gust)
            # Create and solve dynamic problem
            dynamicProblem[i,j,k] = create_DynamicProblem(model=dynamicModel,x0=trimProblem[i,j,k].x[1:end-2],finalTime=tf[i,j],Δt=Δt[i,j],skipInitialStatesUpdate=true,systemSolver=NRdyn)
            solve!(dynamicProblem[i,j,k])
            # Outputs
            t[i,j,k] = dynamicProblem[i,j,k].savedTimeVector
            Nt = length(t[i,j,k])
            rootPitch[i,j,k] = [dynamicProblem[i,j,k].aeroVariablesOverTime[tt][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N for tt in 1:Nt]
            rootAoA[i,j,k] = [dynamicProblem[i,j,k].aeroVariablesOverTime[tt][midSpanElem].flowAnglesAndRates.αₑ-wingAirfoil.attachedFlowParameters.α₀N for tt in 1:Nt]
            Δu2[i,j,k] = [dynamicProblem[i,j,k].nodalStatesOverTime[tt][midSpanElem].u_n2[2] for tt in 1:Nt]
            Δu3[i,j,k] = [dynamicProblem[i,j,k].nodalStatesOverTime[tt][midSpanElem].u_n2[3] for tt in 1:Nt] .- dynamicProblem[i,j,k].nodalStatesOverTime[1][midSpanElem].u_n2[3]
            airspeed[i,j,k] = [dynamicProblem[i,j,k].aeroVariablesOverTime[tt][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for tt in 1:Nt]
            cn[i,j,k] = [dynamicProblem[i,j,k].aeroVariablesOverTime[tt][midSpanElem].aeroCoefficients.cn for tt in 1:Nt]
            cm[i,j,k] = [dynamicProblem[i,j,k].aeroVariablesOverTime[tt][midSpanElem].aeroCoefficients.cm for tt in 1:Nt]
            root_κ1[i,j,k] = [dynamicProblem[i,j,k].compElementalStatesOverTime[tt][midSpanElem].κ[1] for tt in 1:Nt]
            root_κ2[i,j,k] = [dynamicProblem[i,j,k].compElementalStatesOverTime[tt][midSpanElem].κ[2] for tt in 1:Nt]
            tip_Δu3[i,j,k] = [(dynamicProblem[i,j,k].nodalStatesOverTime[tt][1].u_n1[3] - dynamicProblem[i,j,k].nodalStatesOverTime[tt][midSpanElem].u_n2[3]) for tt in 1:Nt]
        end
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosDarpaGustDurationRange"
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
τColors = cgrad(:rainbow, length(τRange), categorical=true)
γstr = round(Int,γ*100)

# Root angle of attack across gust durations
plt_AoA_ate = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,maximum(tf)], ylims=[-15,90], yticks=vcat(-15:15:90), tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plt_AoA_dse = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,maximum(tf)], ylims=[-15,90], yticks=vcat(-15:15:90), tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (k,τ) in enumerate(τRange)
    plot!(plt_AoA_ate,t[1,1,k], rootAoA[1,1,k]*180/π, c=τColors[k], lw=lw, label=string("\$H\$ = ", round(Int,U*τ/(4*.3048))))
    plot!(plt_AoA_dse,t[2,1,k], rootAoA[2,1,k]*180/π, c=τColors[k], lw=lw, label=false)
end
display(plt_AoA_ate)
display(plt_AoA_dse)

# Root bending curvature across gust durations
plt_k2_ate = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", xlims=[0,maximum(tf)], ylims=[-Inf,0], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_k2_dse = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", xlims=[0,maximum(tf)], ylims=[-Inf,0], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (k,τ) in enumerate(τRange)
    plot!(plt_k2_ate,t[1,1,k], root_κ2[1,1,k], c=τColors[k], lw=lw, label=false)
    plot!(plt_k2_dse,t[2,1,k], root_κ2[2,1,k], c=τColors[k], lw=lw, label=false)
end
display(plt_k2_ate)
display(plt_k2_dse)

# Root torsional curvature across gust durations
plt_k1_ate = plot(xlabel="Time [s]", ylabel="Root torsional curvature [1/m]", xlims=[0,maximum(tf)], ylims=[-1e-3,1e-3], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_k1_dse = plot(xlabel="Time [s]", ylabel="Root torsional curvature [1/m]", xlims=[0,maximum(tf)], ylims=[-1e-3,1e-3], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (k,τ) in enumerate(τRange)
    plot!(plt_k1_ate,t[1,1,k], root_κ1[1,1,k], c=τColors[k], lw=lw, label=false)
    plot!(plt_k1_dse,t[2,1,k], root_κ1[2,1,k], c=τColors[k], lw=lw, label=false)
end
display(plt_k1_ate)
display(plt_k1_dse)

# Root bending curvature vs torsional curvature across gust durations
plt_k1k2_ate = plot(xlabel="Root torsional curvature [1/m]", ylabel="Root bending curvature [1/m]", ylims=[-Inf,0], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_k1k2_dse = plot(xlabel="Root torsional curvature [1/m]", ylabel="Root bending curvature [1/m]", ylims=[-Inf,0], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (k,τ) in enumerate(τRange)
    plot!(plt_k1k2_ate, root_κ1[1,1,k], root_κ2[1,1,k], c=τColors[k], lw=lw, label=false)
    plot!(plt_k1k2_dse, root_κ1[2,1,k], root_κ2[2,1,k], c=τColors[k], lw=lw, label=false)
end
display(plt_k1k2_ate)
display(plt_k1k2_dse)

if saveFigures
    τistr = round(τRange[1],digits=1)
    τfstr = round(τRange[end],digits=1)
    savefig(plt_AoA_ate,string(absPath,"/heliosDarpaGustDurationRange_AoAate_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_AoA_dse,string(absPath,"/heliosDarpaGustDurationRange_AoADSe_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_k2_ate,string(absPath,"/heliosDarpaGustDurationRange_k2ate_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_k2_dse,string(absPath,"/heliosDarpaGustDurationRange_k2DSe_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_k1_ate,string(absPath,"/heliosDarpaGustDurationRange_k1ate_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_k1_dse,string(absPath,"/heliosDarpaGustDurationRange_k1DSe_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_k1k2_ate,string(absPath,"/heliosDarpaGustDurationRange_k1k2ate_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
    savefig(plt_k1k2_dse,string(absPath,"/heliosDarpaGustDurationRange_k1k2DSe_P",round(P),"_tau",τistr,"_to_",τfstr,"_gamma",γstr,".pdf"))
end

# Loop gust durations
for (k,τ) in enumerate(τRange)

    gr()

    # Root angle of attack
    plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,maximum(tf)], ylims=[-15,90], yticks=vcat(-15:15:90), tickfont=font(ts), guidefont=font(fs), legend=(0.6,0.9), legendfontsize=lfs)
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], rootAoA[i,j,k]*180/π, c=colors[i], ls=ls[j], lw=lw, label=labels[i,j])
        end
    end
    display(plt_AoA)

    # Airspeed
    plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [ft/s]", xlims=[0,maximum(tf)], ylims=[0,Inf], yticks=vcat(0:20:1e3), tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], airspeed[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_U)

    # Elevation change
    plt_alt = plot(xlabel="Time [s]", ylabel="Elevation change [m]", xlims=[0,maximum(tf)], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], Δu3[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_alt)

    # Root bending curvature
    plt_k2 = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", xlims=[0,maximum(tf)], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], root_κ2[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_k2)

    # Tip OOP displacement
    L = 0.3048*40*(2+cosd(10))
    plt_OOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,maximum(tf)], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], tip_Δu3[i,j,k]/L*100, c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_OOP)

    # cn x time
    plt_cnt = plot(xlabel="Time [s]", ylabel="\$c_n\$", xlims=[0,maximum(tf)], ylims=[-0.5,3], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], cn[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_cnt)

    # cm  x time
    plt_cmt = plot(xlabel="Time [s]", ylabel="\$c_m\$", xlims=[0,maximum(tf)], ylims=[-0.6,0.2], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(t[i,j,k], cm[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_cmt)

    # cn x α
    plt_cna = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_n\$", xlims=[-Inf,Inf], ylims=[-Inf,Inf], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(rootPitch[i,j,k]*180/π, cn[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_cna)

    # cm  x α
    plt_cma = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_m\$", xlims=[-Inf,Inf], ylims=[-Inf,Inf], tickfont=font(ts), guidefont=font(fs))
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,λ) in enumerate(λRange)
            plot!(rootPitch[i,j,k]*180/π, cm[i,j,k], c=colors[i], ls=ls[j], lw=lw, label=false)
        end
    end
    display(plt_cma)

    # Strings
    τstr = round(τ,digits=1)
    γstr = round(Int,γ*100)

    # Save figures, if applicable
    if saveFigures
        savefig(plt_AoA,string(absPath,"/heliosDarpaGustDurationRange_AoA_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_U,string(absPath,"/heliosDarpaGustDurationRange_U_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_alt,string(absPath,"/heliosDarpaGustDurationRange_alt_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_k2,string(absPath,"/heliosDarpaGustDurationRange_k2_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_OOP,string(absPath,"/heliosDarpaGustDurationRange_OOP_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_cnt,string(absPath,"/heliosDarpaGustDurationRange_cnt_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_cmt,string(absPath,"/heliosDarpaGustDurationRange_cmt_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_cna,string(absPath,"/heliosDarpaGustDurationRange_cna_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
        savefig(plt_cma,string(absPath,"/heliosDarpaGustDurationRange_cma_P",round(P),"_tau",τstr,"_gamma",γstr,".pdf"))
    end

    # Maximum altitude and forward distance changes for animation plot limits
    Δforward_11 = round(Int,maximum(Δu2[1,1,k]))
    Δforward_21 = round(Int,maximum(Δu2[2,1,k]))
    Δalt_min_11 = min(-20,round(Int,minimum(Δu3[1,1,k])))
    Δalt_max_11 = max(+20,round(Int,maximum(Δu3[1,1,k])))
    Δalt_min_21 = min(-20,round(Int,minimum(Δu3[2,1,k])))
    Δalt_max_21 = max(+20,round(Int,maximum(Δu3[2,1,k])))

    # Animations
    if createAnimations
        plot_dynamic_deformation(dynamicProblem[1,1,k],refBasis="I",view=(60,15),plotDistLoads=false,plotFrequency=2,plotLimits=([-40,40],[-5,Δforward_11],[Δalt_min_11,Δalt_max_11]),followAssembly=true,save=saveFigures,savePath=string(relPath,"/heliosDarpaGustDurationRange_attached_elastic_P",round(P),"_tau",τstr,"_gamma",γstr,".gif"),displayProgress=true)
        plot_dynamic_deformation(dynamicProblem[2,1,k],refBasis="I",view=(60,15),plotDistLoads=false,plotFrequency=2,plotLimits=([-40,40],[-5,Δforward_21],[Δalt_min_21,Δalt_max_21]),followAssembly=true,save=saveFigures,savePath=string(relPath,"/heliosDarpaGustDurationRange_ds_elastic_P",round(P),"_tau",τstr,"_gamma",γstr,".gif"),displayProgress=true)
    end
end

println("Finished heliosDarpaGustDurationRange.jl")