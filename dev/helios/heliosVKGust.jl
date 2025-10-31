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
seeds = [10]
τ = 360
γ = 0.15
t0 = 1
tAfter = 120
gustResolution = 0.005 # [m]
σ = U*γ
gustLength = U*τ
N = ceil(Int,gustLength/gustResolution+1)
c0 = [0; U*t0; 0]

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
Δu1 = Array{Vector{Float64}}(undef,length(seeds),length(aeroSolvers))
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
relPath = "/dev/helios/figures/heliosVKGust"
absPath = string(pwd(),relPath)
mkpath(absPath)
mkpath(string(pwd(),"/dev/helios/data/heliosVKGust"))

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
fps = 15
labels = ["Attached-flow" "Dynamic stall"]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)

# Sweep seeds
for (i,seed) in enumerate(seeds)
    # Generate gust 
    gust = create_Continuous1DSpaceGust(spectrum="vK",gustLength=gustLength,N=N,σ=σ,c0=c0,seed=seed,plotPSD=false)
    # Sweep aero solvers
    for (j,aeroSolver) in enumerate(aeroSolvers)
        println("Solving for seed $seed, $(aeroSolver.name) solver")
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
        Δu1[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[1] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[1]
        Δu2[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[2] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[2]
        Δu3[i,j] = [dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3] for k in 1:Nt] .- dynamicProblem[i,j].nodalStatesOverTime[1][midSpanElem].u_n2[3]
        airspeed[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].flowVelocitiesAndRates.U∞/.3048 for k in 1:Nt]
        cn[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cn for k in 1:Nt]
        cm[i,j] = [dynamicProblem[i,j].aeroVariablesOverTime[k][midSpanElem].aeroCoefficients.cm for k in 1:Nt]
        root_κ1[i,j] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[1] for k in 1:Nt]
        root_κ2[i,j] = [dynamicProblem[i,j].compElementalStatesOverTime[k][midSpanElem].κ[2] for k in 1:Nt]
        tip_Δu3[i,j] = [(dynamicProblem[i,j].nodalStatesOverTime[k][1].u_n1[3] - dynamicProblem[i,j].nodalStatesOverTime[k][midSpanElem].u_n2[3]) for k in 1:Nt]
        # Find times where gust begins and ends
        tig = round(Int,t0/Δt+1)
        tfg = findfirst(x -> x >= U*(τ+t0), U*t[i,j] .+ Δu2[i,j])
        if isnothing(tfg)
            tfg = Nt
        end
        gustDur[i,j] = t[i,j][tfg] - t[i,j][tig]
        # Aircraft semispan [m]
        L = 0.3048*40*(2+cosd(10))
        # Arrays during gust
        timeInGust = t[i,j][tig:tfg]
        rootAoAInGust = rootAoA[i,j][tig:tfg]
        airspeedInGust = airspeed[i,j][tig:tfg]
        rootκ1InGust = root_κ1[i,j][tig:tfg]
        rootκ2InGust = root_κ2[i,j][tig:tfg]
        tipu3InGust = tip_Δu3[i,j][tig:tfg]
        # Compute counts of exceedance
        COE_rootAoA[i,j] = count_of_exceedance(time=timeInGust, series=(rootAoAInGust.-rootAoAInGust[1])*180/π, datum=0, marker=rootAoAMarkers)
        COE_rootκ1[i,j] = count_of_exceedance(time=timeInGust, series=(rootκ1InGust.-rootκ1InGust[1]), datum=0, marker=root_κ1Markers)
        COE_rootκ2[i,j] = count_of_exceedance(time=timeInGust, series=-(rootκ2InGust.-rootκ2InGust[1]), datum=0, marker=root_κ2Markers)
        COE_tipu3[i,j] = count_of_exceedance(time=timeInGust, series=(tipu3InGust.-tipu3InGust[1])/L*100, datum=0, marker=tip_u3Markers)
        # Save
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("COE_rootAoA_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_rootAoA[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("COE_rootκ1_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_rootκ1[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("COE_rootκ2_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_rootκ2[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("COE_tipu3_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",COE_tipu3[i,j])
        writedlm(pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("gustDur_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt",gustDur[i,j])

        @save pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("timeInGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" timeInGust
        @save pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("rootAoAInGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" rootAoAInGust
        @save pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("airspeedInGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" airspeedInGust
        @save pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("rootκ1InGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" rootκ1InGust
        @save pkgdir(AeroBeams)*"/dev/helios/data/heliosVKGust/"*string("rootκ2InGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" rootκ2InGust
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
    
    # Spanwise-direction change
    plt_deltax = plot(xlabel="Time [s]", ylabel="Spanwise-direction change [m]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], Δu1[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_deltax)

    # Flight-direction change
    plt_deltay = plot(xlabel="Time [s]", ylabel="Forward-direction change [m]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(t[i,j], Δu2[i,j], c=colors[j], lw=lw, label=false)
    end
    display(plt_deltay)

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

    # Plot gust velocity components
    X = c0[2] .+ LinRange(0,gustLength,N)
    Xarray = [[0; x; 0] for x in X]
    u = [gust.U(Xarray[i]) for i in eachindex(Xarray)]
    v = [gust.V(Xarray[i]) for i in eachindex(Xarray)]
    w = [gust.W(Xarray[i]) for i in eachindex(Xarray)]
    plt_uvw = plot(xlabel="\$x_2\$ [m]", ylabel="Gust velocity [m/s]",tickfont=font(ts),guidefont=font(fs),legendfontsize=12,legend=:topright, xlims=[0,4500])
    plot!(X,u, lw=2, c=:purple, label="\$u\$")
    plot!(X,v, lw=2, c=:green, label="\$v\$")
    plot!(X,w, lw=2, c=:red, label="\$w\$")
    display(plt_uvw)
    savefig(plt_uvw,string(absPath,"/heliosVKGust_seed",seed,"_gustVelocities.pdf"))

    # Save figures, if applicable
    if saveFigures
        savefig(plt_AoA,string(absPath,"/heliosVKGust_AoA_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_U,string(absPath,"/heliosVKGust_U_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_deltax,string(absPath,"/heliosVKGust_deltax_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_deltay,string(absPath,"/heliosVKGust_deltay_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_deltaz,string(absPath,"/heliosVKGust_deltaz_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_k1,string(absPath,"/heliosVKGust_k1_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_k2,string(absPath,"/heliosVKGust_k2_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_OOP,string(absPath,"/heliosVKGust_OOP_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cnt,string(absPath,"/heliosVKGust_cnt_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cmt,string(absPath,"/heliosVKGust_cmt_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cna,string(absPath,"/heliosVKGust_cna_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
        savefig(plt_cma,string(absPath,"/heliosVKGust_cma_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".pdf"))
    end

    # Maximum and minimum position changes for animation plot limits
    Δu1_min = min(-40,round(Int,minimum(Δu1[i,1])),round(Int,minimum(Δu1[i,2])))
    Δu1_max = max(40,round(Int,maximum(Δu1[i,1])),round(Int,maximum(Δu1[i,2])))
    Δu2_min = min(-40,round(Int,minimum(Δu2[i,1])),round(Int,minimum(Δu2[i,2])))
    Δu2_max = max(40,round(Int,maximum(Δu2[i,1])),round(Int,maximum(Δu2[i,2])))
    Δu3_min = min(-40,round(Int,minimum(Δu3[i,1])),round(Int,minimum(Δu3[i,2])))
    Δu3_max = max(40,round(Int,maximum(Δu3[i,1])),round(Int,maximum(Δu3[i,2])))

    # Animation
    Δtanim = 1
    anim = plot_dynamic_deformations(dynamicProblem[i,:],refBasis="I",view=(60,15),legendEntries=["Attached-flow", "Dynamic stall"],plotBCs=false,plotDistLoads=false,plotFrequency=ceil(Int,Δtanim/Δt),fps=fps,followAssembly=true,plotLimits=([Δu1_min,Δu1_max],[Δu2_min,Δu2_max],[Δu3_min,Δu3_max]),save=true,savePath=string(relPath,"/heliosVKGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,".gif"),displayProgress=true)
    display(anim)

    # GC
    # dynamicProblem[i,1] = dynamicProblem[i,2] = create_DynamicProblem(model=dummyDynamicModel,finalTime=tf,Δt=Δt)
    # GC.gc()
end

println("Finished heliosVKGust.jl")