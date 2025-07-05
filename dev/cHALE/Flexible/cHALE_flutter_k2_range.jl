using AeroBeams, LinearInterpolations, JLD2, Base.Threads

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 2

# Altitude
h = 20e3

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Discretization
if λ == 1
    nElemWing = 80
elseif λ > 1
    nElemWing = 40
end
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Set number of vibration modes
nModes = 25

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
if λ == 1
    URange = unique(sort(vcat(20:0.5:50,42.3,42.4)))
elseif λ == 2
    URange = unique(sort(vcat(20:0.5:60)))
end

# Damping ratio shift for T-IP modes (in order to match flutter speeds of dynamic solutions)
if λ == 1
    σShift = 3e-2
elseif λ == 2
    σShift = 1.5e-2
end

# Damping ratio tolerance for flutter detection at lower speed limit
σtol = 2e-3

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
trimProblemSpringed = Array{TrimProblem}(undef,length(k2Range),length(URange))
eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))

trimAoA = fill(NaN, length(k2Range), length(URange))
trimThrust = fill(NaN, length(k2Range), length(URange))
trimδ = fill(NaN, length(k2Range), length(URange))
trimEnduranceFactor = fill(NaN, length(k2Range), length(URange))

untrackedFreqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
damps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeDampings = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]
modeDampingRatios = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]

x1_0 = Array{Vector{Float64}}(undef,length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(k2Range))
x1_e_wing = Array{Vector{Float64}}(undef,length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
twist = Array{Vector{Float64}}(undef,length(k2Range),length(URange))

# Attachment springs' stiffness
μu = 1e-2
μp = 1e-2

# ELement ranges
elemRangeRightWing = 1 + div(nElemWing,2) : nElemWing

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Set attachment springs
    spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])
    spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        # Create and trim problem
        trimProblem[i,j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        trimThrust[i,j] = stabilizersAero ? trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling : trimProblem[i,j].x[end]*trimProblem[i,j].model.forceScaling
        trimδ[i,j] = stabilizersAero ? trimProblem[i,j].x[end] : 0
        println("Trim AoA = $(trimAoA[i,j]*180/π), trim thrust = $(trimThrust[i,j]), trim δ = $(trimδ[i,j]*180/π)")
        lift = trimProblem[i,j].model.mass*trimProblem[i,j].model.atmosphere.g - trimThrust[i,j]*sin(trimAoA[i,j])
        drag = trimThrust[i,j]*cos(trimAoA[i,j])
        qS = 1/2*trimProblem[i,j].model.atmosphere.ρ*U^2*(32*1)
        trimEnduranceFactor[i,j] = lift^1.5/drag*sqrt(1/qS)
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(cHALEtrim.elements[e].r_n1[1],cHALEtrim.elements[e].r_n2[1]) for e in elemRangeRightWing]...)
            x3_0[i] = vcat([vcat(cHALEtrim.elements[e].r_n1[3],cHALEtrim.elements[e].r_n2[3]) for e in elemRangeRightWing]...)
            # Undeformed elemental positions
            x1_e_wing[i] = [cHALEtrim.elements[e].x1 for e in elemRangeRightWing]
        end
        # Displacements over span
        u1_of_x1[i,j] = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in elemRangeRightWing]...)
        u3_of_x1[i,j] = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in elemRangeRightWing]...)
        u1_of_x1[i,j] .-= u1_of_x1[i,j][1]
        u3_of_x1[i,j] .-= u3_of_x1[i,j][1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] .+ u1_of_x1[i,j]
        x3_def[i,j] = x3_0[i] .+ u3_of_x1[i,j]
        # Angle of twist
        p1_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[1]) for e in elemRangeRightWing]...)
        p2_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[2],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[2]) for e in elemRangeRightWing]...)
        p3_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[3]) for e in elemRangeRightWing]...)
        twist[i,j] = [asind((first(rotation_tensor_WM([p1_of_x1[k],p2_of_x1[k],p3_of_x1[k]]))*AeroBeams.a2)[3]) for k in eachindex(p1_of_x1)]
        twist[i,j] .-= twist[i,j][1] # discount root angle (rigid-body rotation)
        # Model for trim problem with springs
        cHALEtrimSpringed,_,_,tailBoomSpringed,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
        # Add springs
        add_springs_to_beam!(beam=tailBoomSpringed,springs=[spring1,spring2])
        # Update model
        cHALEtrimSpringed.skipValidationMotionBasisA = true
        update_model!(cHALEtrimSpringed)
        # Create and solve trim problem with springs
        trimProblemSpringed[i,j] = create_TrimProblem(model=cHALEtrimSpringed,systemSolver=NRtrim,x0=trimProblem[i,j].x)
        solve!(trimProblemSpringed[i,j])
        # Retrieve and compare trim outputs
        trimAoASpringed = (trimProblemSpringed[i,j].aeroVariablesOverσ[end][cHALEtrimSpringed.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblemSpringed[i,j].aeroVariablesOverσ[end][cHALEtrimSpringed.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
        trimThrustSpringed = trimProblemSpringed[i,j].x[end-1]*trimProblemSpringed[i,j].model.forceScaling
        trimδSpringed = trimProblemSpringed[i,j].x[end]
        println("Trim outputs ratios springed/nominal: AoA = $(trimAoASpringed/trimAoA[i,j]), T = $(trimThrustSpringed/trimThrust[i,j]), δ = $(trimδSpringed/trimδ[i,j])")
        # Model for eigen problem
        cHALEeigen,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδSpringed,thrust=trimThrustSpringed,k2=k2,hasInducedDrag=hasInducedDrag)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=cHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf],refTrimProblem=trimProblemSpringed[i,j])
        solve_eigen!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Frequencies and dampings after mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

# Adjusted dampings with damping ratio shift for T-IP modes
if λ == 1
    TIP_modes = [[9,14], [9,15], [10,13], [12,13], [13,11]]
elseif λ == 2
    TIP_modes = [[9,14], [11,14], [11,14], [11,15], [10,15]]
end
modeDampingsAdj = deepcopy(modeDampings)
modeDampingRatiosAdj = deepcopy(modeDampingRatios)
for (i,k2) in enumerate(k2Range)
    for mode in TIP_modes[i]
        modeDampingRatiosAdj[i,mode] .-= σShift
        modeDampingsAdj[i,mode] .-= σShift*modeFrequencies[i,mode]
    end
end

# Compute flutter variables
if λ == 1
    modes2ignoreFlutter = [[], [], [], [4], [4]] # Dutch-roll mode
elseif λ == 2
    modes2ignoreFlutter = [[], [], [], [], []] # Dutch-roll mode
end
flutterOnsetMode = [fill(0,0) for _ in 1:length(k2Range)]
flutterOffsetMode = [fill(0,0) for _ in 1:length(k2Range)]
flutterOnsetSpeedOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOnsetFreqOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOnsetTipOOPOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOffsetSpeedOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOffsetFreqOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOnsetSpeed = [Inf for _ in 1:length(k2Range)]
flutterOnsetFreq = [NaN for _ in 1:length(k2Range)]
flutterOnsetTipOOP = [NaN for _ in 1:length(k2Range)]
flutterOnsetSpeedsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
flutterOnsetFreqsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
flutterOffsetSpeedsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
flutterOffsetFreqsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
indicesFlutterOnset = [fill(0,0) for _ in 1:length(k2Range)]
for (i,k2) in enumerate(k2Range)
    # Flutter onset/offset data of each mode
    for mode in 1:nModes
        if mode in modes2ignoreFlutter[i]
            continue
        end
        iOnset = findall(j -> modeDampingRatiosAdj[i,mode][j] < 0 && modeDampingRatiosAdj[i,mode][j+1] > 0, 1:length(URange)-1)
        iOffset = findall(j -> modeDampingRatiosAdj[i,mode][j] > 0 && modeDampingRatiosAdj[i,mode][j+1] < 0, 1:length(URange)-1)
        if modeDampingRatiosAdj[i,mode][1] > σtol
            push!(flutterOnsetMode[i],mode)
            push!(flutterOnsetSpeedOfMode[i,mode],URange[1])
            push!(flutterOnsetFreqOfMode[i,mode],modeFrequencies[i,mode][1])
            push!(flutterOnsetTipOOPOfMode[i,mode],x3_def[i,1][end]/16*100)
        end
        if !isempty(iOnset)
            for iO in iOnset
                push!(flutterOnsetMode[i],mode)
                push!(flutterOnsetSpeedOfMode[i,mode],interpolate(modeDampingRatiosAdj[i,mode][iO:iO+1],URange[iO:iO+1],0))
                push!(flutterOnsetFreqOfMode[i,mode],interpolate(modeDampingRatiosAdj[i,mode][iO:iO+1],modeFrequencies[i,mode][iO:iO+1],0))
                push!(flutterOnsetTipOOPOfMode[i,mode],(x3_def[i,iO][end]+x3_def[i,iO+1][end])/2/16*100)
            end
        end
        if !(isempty(iOffset) || isempty(iOnset))
            for iO in iOffset
                push!(flutterOffsetMode[i],mode)
                push!(flutterOffsetSpeedOfMode[i,mode],interpolate(-modeDampingRatiosAdj[i,mode][iO:iO+1],URange[iO:iO+1],0))
                push!(flutterOffsetFreqOfMode[i,mode],interpolate(-modeDampingRatiosAdj[i,mode][iO:iO+1],modeFrequencies[i,mode][iO:iO+1],0))
            end
        end
    end
    # All flutter onset/offset speeds and frequencies and tip OOP
    flutterOnsetSpeedsAll[i] = vcat(filter(!isempty,flutterOnsetSpeedOfMode[i,:])...)
    flutterOnsetFreqsAll[i] = vcat(filter(!isempty,flutterOnsetFreqOfMode[i,:])...)
    flutterOffsetSpeedsAll[i] = vcat(filter(!isempty,flutterOffsetSpeedOfMode[i,:])...)
    flutterOffsetFreqsAll[i] = vcat(filter(!isempty,flutterOffsetFreqOfMode[i,:])...)
    flutterOnsetTipOOPAll = vcat(filter(!isempty,flutterOnsetTipOOPOfMode[i,:])...)
    # Ordered indices (of URange) of flutter onset
    flutterOnsetSpeedsAllOrdered = sort(flutterOnsetSpeedsAll[i])
    for speed in flutterOnsetSpeedsAllOrdered
        push!(indicesFlutterOnset[i], -1 + findfirst(x-> x > speed, URange))
    end
    # Lowest flutter onset speed, corresponding frequency and tip OOP
    if !isempty(flutterOnsetSpeedsAll[i])
        iLowest = sortperm(flutterOnsetSpeedsAll[i])[1]
        flutterOnsetSpeed[i] = flutterOnsetSpeedsAll[i][iLowest]
        flutterOnsetFreq[i] = flutterOnsetFreqsAll[i][iLowest]
        flutterOnsetTipOOP[i] = flutterOnsetTipOOPAll[iLowest]
    end
end

# Set paths
relPathFig = "/dev/cHALE/Flexible/outputs/figures/cHALE_flutter_k2_range"
absPathFig = string(pwd(),relPathFig)
mkpath(absPathFig)

relPathData = "/dev/cHALE/Flexible/outputs/data/cHALE_flutter_k2_range/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)

absPathDataWingMatched = string(pwd(),"/dev/cHALE/Flexible/outputs/data/cHALEwing_flutter_matchedPitch_k2_range/")
absPathDataWing = string(pwd(),"/dev/cHALE/Flexible/outputs/data/cHALEwing_g0_flutter_fixedPitch_k2_range/")

# Save and load flutter data
flutterOnsetSpeedAircraft = flutterOnsetSpeed
flutterOnsetFreqAircraft = flutterOnsetFreq
flutterOnsetTipOOPAircraft = flutterOnsetTipOOP
@save absPathData*string("aircraft_lambda",λ,"_allFlutterSpeed.jld2") flutterOnsetSpeedsAll
@save absPathData*string("aircraft_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedAircraft
@save absPathData*string("aircraft_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqAircraft
@save absPathData*string("aircraft_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPAircraft

@load absPathDataWingMatched*string("wing_matched_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedWingMatched
@load absPathDataWingMatched*string("wing_matched_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqWingMatched
@load absPathDataWingMatched*string("wing_matched_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPWingMatched
@load absPathDataWing*string("wing_g0_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedWingg0
@load absPathDataWing*string("wing_g0_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqWingg0
@load absPathDataWing*string("wing_g0_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPWingg0

using Plots, ColorSchemes

# # Mode shapes at lowest airspeed
# for (n,k2) in enumerate(k2Range)
#     plt = plot_mode_shapes(eigenProblem[n,1],nModes=6,scale=5,view=(30,30),modalColorScheme=:rainbow,legendPos=:outertop,save=true,savePath=string(relPathFig,string("/cHALE_flutter_k2_range_modeShapes_k2",k2,"_lambda",λ,".pdf")))
#     display(plt)
# end

# Plot configurations
colors = cgrad(:rainbow, length(k2Range), categorical=true)
ts = 10
fs = 16
lfs = 10
tsz = 10
lw = 2
ms = 3
msw = 0
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,250], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampingsAdj[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampingsAdj[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_rootlocus_lambda",λ,".pdf"))

# Root locus (zoom on low frequency modes)
if λ == 1
    dampLim = [-10,2]
    freqLim = [0,50]
elseif λ == 2
    dampLim = [-10,2]
    freqLim = [0,80]
end
plt_RLlow = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=(0.12,0.8))
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampingsAdj[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampingsAdj[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
if λ == 1
    phugoidTextPos = [1.25, 1.5]
    DRTextPos = [0.8, 5]
    TIP1textPos = [0.75, 8]
    TIP2textPos = [0.75, 32]
    annotate!(phugoidTextPos[1], phugoidTextPos[2], text("phugoid", tsz))
    annotate!(DRTextPos[1], DRTextPos[2], text("D-R", tsz))
    annotate!(TIP1textPos[1], TIP1textPos[2], text("1st T-IP", tsz))
    annotate!(TIP2textPos[1], TIP2textPos[2], text("2nd T-IP", tsz))
    quiver!([phugoidTextPos[1]-0.7,DRTextPos[1]-0.35,TIP1textPos[1]-0.25,TIP2textPos[1]-0.25], [phugoidTextPos[2]-0.75,DRTextPos[2]-0.7,TIP1textPos[2]+1.5,TIP2textPos[2]-1.5], quiver=([-0.35,-0.35,-0.25,-0.35], [-0.3,-1.8,3,-2.2]), arrow=:closed, linecolor=:black)
end
display(plt_RLlow)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_rootlocus_low_lambda",λ,".pdf"))

# Root locus (zoom on T-IP modes)
if λ == 1
    dampLim = [-2.5,2]
    freqLim = [2,45]
elseif λ == 2
    dampLim = [-4,4]
    freqLim = [5,80]
end
plt_RLTIP = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    for mode in TIP_modes[i]
        scatter!(modeDampingsAdj[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampingsAdj[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
if λ == 1
    TIP1textPos = [0.5, 5]
    TIP2textPos = [0.5, 35]
    annotate!(TIP1textPos[1], TIP1textPos[2], text("1st T-IP", 12))
    annotate!(TIP2textPos[1], TIP2textPos[2], text("2nd T-IP", 12))
    quiver!([TIP1textPos[1]-0.35,TIP2textPos[1]-0.35], [TIP1textPos[2],TIP2textPos[2]], quiver=([-0.65,-0.55], [3,-2]), arrow=:closed, linecolor=:black)
end
display(plt_RLTIP)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_rootlocus_TIP_lambda",λ,".pdf"))

# Root locus (zoom on phugoid mode)
if λ == 1
    dampLim = [-0.15,0.15]
    freqLim = [0,0.5]
    iUtransition_k3 = findfirst(x->x>34,URange)
    iUtransition_k4 = findfirst(x->x>39,URange)
    iUtransition_k5 = findfirst(x->x>42.3,URange)
    phugoidMode = [fill(2,length(URange)), 
                   fill(2,length(URange)), 
                   vcat(fill(1,1:iUtransition_k3),fill(nModes,length(URange)-iUtransition_k3)), 
                   vcat(fill(2,1:iUtransition_k4),fill(nModes,length(URange)-iUtransition_k4)), 
                   vcat(fill(2,1:iUtransition_k5),fill(nModes,length(URange)-iUtransition_k5))]
elseif λ == 2
    dampLim = [-0.05,0.2]
    freqLim = [0,0.5]
    phugoidMode = [fill(1,length(URange)), fill(2,length(URange)), fill(2,length(URange)), fill(1,length(URange)), fill(1,length(URange))]
end
plt_RLphugoid = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    for j in eachindex(URange)
        scatter!([modeDampingsAdj[i,phugoidMode[i][j]][j]], [modeFrequencies[i,phugoidMode[i][j]][j]], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    scatter!([modeDampingsAdj[i,phugoidMode[i][1]][1]], [modeFrequencies[i,phugoidMode[i][1]][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RLphugoid)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_rootlocus_phugoid_lambda",λ,".pdf"))

# V-g-f
if λ == 1
    dampLim = [-0.3,0.2]
    freqLim = [0,50]
elseif λ == 2
    dampLim = [-0.3,0.25]
    freqLim = [0,80]
end
for (i,k2) in enumerate(k2Range)
    plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=freqLim, tickfont=font(ts), guidefont=font(12))
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=dampLim, tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
    plot!(plt_Vg,URange,zeros(length(URange)), c=:gray, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URange, modeDampingRatiosAdj[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(string(absPathFig,"/cHALE_flutter_k2_range_Vgf",i,"_lambda",λ,".pdf"))
end

# Flutter onset speeds vs k2
if λ == 1
    ULim = [0,50]
elseif λ == 2
    ULim = [0,65]
end
plt_Uf = plot(xlabel="\$k_2\$ [1/m]", ylabel="Flutter speed [m/s]", ylims=ULim, xticks=k2Range, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomleft)
plot!(k2Range,flutterOnsetSpeedWingg0, marker=(:circle, 10), msw=msw, c=:black, lw=lw, ls=:dash, label="Wing - zero pitch, no gravity")
plot!(k2Range,flutterOnsetSpeedWingMatched, marker=(:square, 10), msw=msw, c=:black, lw=lw, ls=:dot, label="Wing - trim-matched pitch")
plot!(k2Range,flutterOnsetSpeedAircraft, marker=(:diamond, 10), msw=0, c=:black, lw=lw, ls=:solid, label="Aircraft - trimmed flight")
display(plt_Uf)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_speedOn_lambda",λ,".pdf"))

# Flutter onset frequencies vs k2
if λ == 1
    ωLim = [0,30]
elseif λ == 2
    ωLim = [0,40]
end
plt_ωf = plot(xlabel="\$k_2\$ [1/m]", ylabel="Flutter frequency [rad/s]", ylims=ωLim, xticks=k2Range, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomleft)
plot!(k2Range,flutterOnsetFreqWingg0, marker=(:circle, 10), msw=msw, c=:black, lw=lw, ls=:dash, label=false)
plot!(k2Range,flutterOnsetFreqWingMatched, marker=(:square, 10), msw=msw, c=:black, lw=lw, ls=:dot, label=false)
plot!(k2Range,flutterOnsetFreqAircraft, marker=(:diamond, 10), msw=0, c=:black, lw=lw, ls=:solid, label=false)
display(plt_ωf)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_freqOn_lambda",λ,".pdf"))

# Normalized deformed span at lowest and flutter airspeeds
plt_u3 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[-1,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
plot!([NaN], [NaN], ls=:dash, c=:black, lw=lw, label=string("At flutter"))
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/L, x3_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label=false)
    plot!(-x1_def[i,1]/L, x3_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label=false)
    plot!(x1_def[i,indicesFlutterOnset[i][1]]/L, x3_def[i,indicesFlutterOnset[i][1]]/L, ls=:dash, c=colors[i], lw=lw, label=false)
    plot!(-x1_def[i,indicesFlutterOnset[i][1]]/L, x3_def[i,indicesFlutterOnset[i][1]]/L, ls=:dash, c=colors[i], lw=lw, label=false)
end
display(plt_u3)
savefig(string(absPathFig,"/cHALE_flutter_k2_range_disp_lambda",λ,".pdf"))

# Normalized deformed span at lowest and highest airspeeds
plt_u32 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=(0.4,0.9))
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
plot!([NaN], [NaN], ls=:dash, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[end],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/L, x3_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label="\$k_2 = $(k2) \$")
    plot!(x1_def[i,end]/L, x3_def[i,end]/L, ls=:dash, c=colors[i], lw=lw, label=false)
end
display(plt_u32)
savefig(string(absPathFig,"/cHALE_trim_k2_range_disp_lambda",λ,".pdf"))

# Trim root angle of attack
if λ == 1
    AoALim = [-5,20]
elseif λ == 2
    AoALim = [-5,15]
end
plt_trimAoA = plot(xlabel="Airspeed [m/s]", ylabel="Trim root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=AoALim, tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimAoA[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimAoA)
savefig(string(absPathFig,"/cHALE_trim_k2_range_trimAoA_lambda",λ,".pdf"))

# Trim thrust
if λ == 1
    TLim = [0,100]
elseif λ == 2
    TLim = [0,80]
end
plt_trimThrust = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]], ylims=TLim, tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_trimThrust)
savefig(string(absPathFig,"/cHALE_trim_k2_range_trimThrust_lambda",λ,".pdf"))

# Trim elevator deflection
if λ == 1
    δLim = [-50,10]
elseif λ == 2
    δLim = [-30,10]
end
plt_trimDelta = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[URange[1],URange[end]], ylims=δLim, tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimδ[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimDelta)
savefig(string(absPathFig,"/cHALE_trim_k2_range_trimDelta_lambda",λ,".pdf"))

# Trim endurance factor (cl^1.5/cd)
if λ == 1
    eLim = [0,40]
elseif λ == 2
    eLim = [0,40]
end
plt_enduranceFactor = plot(xlabel="Airspeed [m/s]", ylabel="Trim \$c_L^{1.5}/c_D\$", xlims=[URange[1],URange[end]], ylims=eLim, tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimEnduranceFactor[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_enduranceFactor)
savefig(string(absPathFig,"/cHALE_trim_k2_range_enduranceFactor_lambda",λ,".pdf"))

println("Finished cHALE_flutter_k2_range.jl")