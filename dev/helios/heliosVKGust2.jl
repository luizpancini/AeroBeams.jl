using AeroBeams, DelimitedFiles, JLD2

# Flag to save figures
saveFigures = true

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
P = 150

# Airspeed [m/s]
U = 40*0.3048

# Gust variables
seeds = 1:1:1
τ = 360
γ = 0.1
t0 = 1
tAfter = 180
σ = U*γ

# Time variables
Δt = 2.5e-2
tf = t0+τ+tAfter

# Markers for frequencies of exceedance
rootAoAMarkers = vcat(0:1:90)
root_κ1Markers = vcat(0:5e-5:1e-2)
root_κ2Markers = vcat(0:1e-4:5e-2)
tip_u3Markers = vcat(0:1:50)

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

# Create dummy dynamic model
dummyDynamicModel,_ = create_Helios(aeroSolver=aeroSolvers[1],reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,stiffnessFactor=λ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δ=0,thrust=0)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(seeds),length(aeroSolvers))
dynamicProblem = Array{DynamicProblem}(undef,length(seeds),length(aeroSolvers))
t = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
rootPitch = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
rootAoA = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
Δu2 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
Δu3 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
airspeed = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
cn = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
cm = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
root_κ1 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
root_κ2 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
tip_Δu3 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
COE_rootAoA = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
COE_rootκ1 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
COE_rootκ2 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
COE_tipu3 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
gustDur  = Array{Float64}(undef,length(seeds),length(aeroSolvers))

# Strings
Pstr = round(P)
τstr = round(τ,digits=1)
γstr = round(Int,γ*100)

# Set paths
relPath = "/dev/helios/figures/heliosVKGust2"
absPath = string(pwd(),relPath)
mkpath(absPath)
mkpath(string(pwd(),"/dev/helios/data/heliosVKGust2"))

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
labels = ["Attached flow" "Dynamic stall"]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)

# Sweep seeds
for (i,seed) in enumerate(seeds)
    # Generate gust 
    gust = create_Continuous1DGust(spectrum="vK",generationMethod="sinusoids",initialTime=t0,duration=τ,components=[1,2,3],ωmax=20π,Uref=U,σ=σ,seed=seed,plotPSD=false)
    # Sweep aero solvers
    for (j,aeroSolver) in enumerate(aeroSolvers)
        println("Solving for seed $seed, aeroSolver $j")
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
        dynamicProblem[i,j] = create_DynamicProblem(model=dynamicModel,x0=trimProblem[i,j].x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NRdyn)
        solve!(dynamicProblem[i,j])
        # Outputs
        t[i,j] = dynamicProblem[i,j].savedTimeVector
        Nt = length(t[i,j])
        rootPitch[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
        rootAoA[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowAnglesAndRates.αₑ-wingAirfoil.attachedFlowParameters.α₀N for k in 1:Nt]
        Δu2[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[2] for k in 1:Nt]
        Δu3[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[3]
        airspeed[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for k in 1:Nt]
        cn[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cn for k in 1:Nt]
        cm[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cm for k in 1:Nt]
        root_κ1[i,j] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[1] for k in 1:Nt]
        root_κ2[i,j] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[2] for k in 1:Nt]
        tip_Δu3[i,j] = [(dynamicProblem[i,j].nodalStatesOverTime[k][1].u_n1[3] - dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3]) for k in 1:Nt]
        # Gust duration
        gustDur[i,j] = τ
        # Aircraft semispan [m]
        L = 0.3048*40*(2+cosd(10))
        # Compute counts of exceedance
        COE_rootAoA[i,j] = count_of_exceedance(time=t[i,j][tig:tfg], series=(rootAoA[i,j][tig:tfg].-rootAoA[i,j][1])*180/π, datum=0, marker=rootAoAMarkers)
        COE_rootκ1[i,j] = count_of_exceedance(time=t[i,j][tig:tfg], series=(root_κ1[i,j][tig:tfg].-root_κ1[i,j][1]), datum=0, marker=root_κ1Markers)
        COE_rootκ2[i,j] = count_of_exceedance(time=t[i,j][tig:tfg], series=-(root_κ2[i,j][tig:tfg].-root_κ2[i,j][1]), datum=0, marker=root_κ2Markers)
        COE_tipu3[i,j] = count_of_exceedance(time=t[i,j][tig:tfg], series=(tip_Δu3[i,j][tig:tfg].-tip_Δu3[i,j][1])/L*100, datum=0, marker=tip_u3Markers)
        # Save
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust2/"*string("COE_rootAoA_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_rootAoA[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust2/"*string("COE_rootκ1_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_rootκ1[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust2/"*string("COE_rootκ2_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_rootκ2[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust2/"*string("COE_tipu3_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_tipu3[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust2/"*string("gustDur_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",gustDur[i,j])
    end

    # Plots for current seed
    # --------------------------------------------------------------------------
    gr()
    # Root angle of attack
    plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=(0.6,0.9), legendfontsize=lfs)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], rootAoA[i,j]*180/π, c=colors[j], lw=lw, label=labels[j])
    end
    display(plt_AoA)

    # Airspeed
    plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [ft/s]", xlims=[0,tf], ylims=[0,Inf], yticks=vcat(0:20:1e3), tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], airspeed[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_U)

    # Elevation change
    plt_deltaz = plot(xlabel="Time [s]", ylabel="Elevation change [m]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], Δu3[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_deltaz)

    # Root torsional curvature
    plt_k1 = plot(xlabel="Time [s]", ylabel="Root torsional curvature [1/m]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], root_κ1[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_k1)

    # Root bending curvature
    plt_k2 = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], root_κ2[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_k2)

    # Tip OOP displacement
    L = 0.3048*40*(2+cosd(10))
    plt_OOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], tip_Δu3[i,j]/L*100, c=colors[j], lw=lw, label=false)
    end
    display(plt_OOP)

    # cn x time
    plt_cnt = plot(xlabel="Time [s]", ylabel="\$c_n\$", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], cn[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_cnt)

    # cm  x time
    plt_cmt = plot(xlabel="Time [s]", ylabel="\$c_m\$", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], cm[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_cmt)

    # cn x α
    plt_cna = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_n\$", xlims=[-Inf,Inf], ylims=[-Inf,Inf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(rootPitch[i,j]*180/π, cn[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_cna)

    # cm  x α
    plt_cma = plot(xlabel="Steady angle of attack [deg]", ylabel="\$c_m\$", xlims=[-Inf,Inf], ylims=[-Inf,Inf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(rootPitch[i,j]*180/π, cm[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_cma)

    # Save figures, if applicable
    if saveFigures
        savefig(plt_AoA,string(absPath,"/heliosVKGust_AoA_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_U,string(absPath,"/heliosVKGust_U_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_deltaz,string(absPath,"/heliosVKGust_tauz_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_k1,string(absPath,"/heliosVKGust_k1_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_k2,string(absPath,"/heliosVKGust_k2_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_OOP,string(absPath,"/heliosVKGust_OOP_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cnt,string(absPath,"/heliosVKGust_cnt_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cmt,string(absPath,"/heliosVKGust_cmt_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cna,string(absPath,"/heliosVKGust_cna_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cma,string(absPath,"/heliosVKGust_cma_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        # Maximum altitude and forward distance changes for animation plot limits
        Δforward_DS = max(100,round(Int,maximum(Δu2[i,2])))
        Δalt_min_DS = min(-20,round(Int,minimum(Δu3[i,2])))
        Δalt_max_DS = max(+20,round(Int,maximum(Δu3[i,2])))
        # Animations
        # plot_dynamic_deformation(dynamicProblem[i,2],refBasis="I",view=(60,15),plotDistLoads=false,plotFrequency=5,followAssembly=true,plotLimits=([-40,40],[-10,Δforward_DS],[Δalt_min_DS,Δalt_max_DS]),save=saveFigures,savePath=string(relPath,"/heliosVKGust_ds_elastic_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".gif"),displayProgress=true)
    end

    # GC
    dynamicProblem[i,1] = dynamicProblem[i,2] = create_DynamicProblem(model=dummyDynamicModel,finalTime=tf,Δt=Δt)
    GC.gc()
end

println("Finished heliosVKGust2.jl")