using AeroBeams, LinearInterpolations, JLD2

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Stiffness factor
λ = 2

# Altitude
h = 20e3

# Discretization
if λ == 1
    nElemWing = 80
elseif λ > 1
    nElemWing = 40
end
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=1,displayStatus=false)

# Set number of vibration modes
nModes = 5

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
if λ == 1
    URange = unique(sort(vcat(20:0.5:50,42.3,42.4)))
elseif λ == 2
    URange = unique(sort(vcat(20:0.5:65,60:0.25:65)))
end

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))

trimAoA = fill(NaN, length(k2Range), length(URange))

untrackedFreqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
damps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeDampings = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]

highestConvUindex = Array{Int64}(undef,length(k2Range))

x1_0 = Array{Vector{Float64}}(undef, length(k2Range))
x3_0 = Array{Vector{Float64}}(undef, length(k2Range))
x1_n = Array{Vector{Float64}}(undef, length(k2Range))
x1_e = Array{Vector{Float64}}(undef, length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        if λ == 2 && k2 == 0.015 && U ≈ 68.5
            jU325 = findfirst(x->x≈32.5,URange)
            x0Trim = trimProblem[i,jU325].x
        end
        # Create and trim problem
        trimProblem[i,j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Skip if unconverged
        if !trimProblem[i,j].systemSolver.convergedFinalSolution
            highestConvUindex[i] = j-1
            break
        else
            highestConvUindex[i] = j
        end
        # Extract trim variables
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        println("Trim AoA = $(trimAoA[i,j]*180/π)")
        # Model for eigen problem
        wingModel,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA[i,j])
        # Set initial guess solution as previous known solution (use better approximation based on AoA for the discontinuous case of k2=0.045 @U≈42.4)
        x0Eig = j == 1 ? zeros(0) : eigenProblem[i,j-1].x
        if λ == 1 && k2 == 0.045 && U ≈ 42.4
            jU22 = findfirst(x->x≈22,URange)
            x0Eig = eigenProblem[i,jU22].x
        elseif λ == 2 && k2 == 0.015 && U ≈ 63.5
            jU325 = findfirst(x->x≈32.5,URange)
            x0Eig = eigenProblem[i,jU325].x
        end
        # Lower frequency limit
        if λ == 1
            fLower = (k2 == 0.015 && U < 25) ? 1e-3 : 1
        elseif λ == 2
            fLower = 1e-1
        end
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[fLower,Inf],systemSolver=NReigen,x0=x0Eig)
        solve!(eigenProblem[i,j])
        # Skip if unconverged
        if !eigenProblem[i,j].systemSolver.convergedFinalSolution
            highestConvUindex[i] = j-1
            break
        else
            highestConvUindex[i] = j
        end
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(wingModel.elements[e].r_n1[1],wingModel.elements[e].r_n2[1]) for e in 1:div(nElemWing,2)]...)
            x3_0[i] = vcat([vcat(wingModel.elements[e].r_n1[3],wingModel.elements[e].r_n2[3]) for e in 1:div(nElemWing,2)]...)
            # Nodal and elemental arclength positions
            x1_n[i] = vcat([vcat(wingModel.elements[e].x1_n1,wingModel.elements[e].x1_n2) for e in 1:div(nElemWing,2)]...)
            x1_e[i] = [wingModel.elements[e].x1 for e in 1:div(nElemWing,2)]
        end
        # Displacements over span
        u1_of_x1[i,j] = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],eigenProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:div(nElemWing,2)]...)
        u3_of_x1[i,j] = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],eigenProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:div(nElemWing,2)]...)
        u1_of_x1[i,j] .-= u1_of_x1[i,j][1]
        u3_of_x1[i,j] .-= u3_of_x1[i,j][1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] .+ u1_of_x1[i,j]
        x3_def[i,j] = x3_0[i] .+ u3_of_x1[i,j]
    end
    # Frequencies and dampings after mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
    end
end

# Fix eigenvalues manually
if λ == 1
    iU21 = findfirst(x->x≈21, URange)
    modeFrequencies[3,1][iU21] = 0
    modeDampings[3,1][iU21] = (modeDampings[3,1][iU21-1]+modeDampings[3,1][iU21+1])/2
elseif λ == 2
    iU245 = findfirst(x->x≈24.5, URange)
    modeFrequencies[1,1][iU245] = 0
    modeDampings[1,1][iU245] = (modeDampings[1,1][iU245-1]+modeDampings[1,1][iU245+1])/2
end

# Compute flutter variables
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
        iOnset = findall(j -> modeDampings[i,mode][j] < 0 && modeDampings[i,mode][j+1] > 0, 1:length(URange)-1)
        iOffset = findall(j -> modeDampings[i,mode][j] > 0 && modeDampings[i,mode][j+1] < 0, 1:length(URange)-1)
        if modeDampings[i,mode][1]/modeFrequencies[i,mode][1] > 0
            push!(flutterOnsetMode[i],mode)
            push!(flutterOnsetSpeedOfMode[i,mode],URange[1])
            push!(flutterOnsetFreqOfMode[i,mode],modeFrequencies[i,mode][1])
            push!(flutterOnsetTipOOPOfMode[i,mode],x3_def[i,1][end]/16*100)
        end
        if !isempty(iOnset)
            for iO in iOnset
                push!(flutterOnsetMode[i],mode)
                push!(flutterOnsetSpeedOfMode[i,mode],interpolate(modeDampings[i,mode][iO:iO+1],URange[iO:iO+1],0))
                push!(flutterOnsetFreqOfMode[i,mode],interpolate(modeDampings[i,mode][iO:iO+1],modeFrequencies[i,mode][iO:iO+1],0))
                push!(flutterOnsetTipOOPOfMode[i,mode],(x3_def[i,iO][end]+x3_def[i,iO+1][end])/2/16*100)
            end
        end
        if !(isempty(iOffset) || isempty(iOnset))
            for iO in iOffset
                push!(flutterOffsetMode[i],mode)
                push!(flutterOffsetSpeedOfMode[i,mode],interpolate(-modeDampings[i,mode][iO:iO+1],URange[iO:iO+1],0))
                push!(flutterOffsetFreqOfMode[i,mode],interpolate(-modeDampings[i,mode][iO:iO+1],modeFrequencies[i,mode][iO:iO+1],0))
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
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_flutter_matchedPitch_k2_range.jl"
absPath = string(pwd(),relPath)
mkpath(absPath)
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALEwing_flutter_matchedPitch_k2_range/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)
absPathDataWing = string(pwd(),"/dev/cHALE/Flexible/outputs/data/cHALEwing_g0_flutter_fixedPitch_k2_range/")

# Save and load flutter data
flutterOnsetSpeedWingMatched = flutterOnsetSpeed
flutterOnsetFreqWingMatched = flutterOnsetFreq
flutterOnsetTipOOPWingMatched = flutterOnsetTipOOP
@save absPathData*string("wing_matched_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedWingMatched
@save absPathData*string("wing_matched_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqWingMatched
@save absPathData*string("wing_matched_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPWingMatched
@load absPathDataWing*string("wing_g0_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedWingg0
@load absPathDataWing*string("wing_g0_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqWingg0
@load absPathDataWing*string("wing_g0_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPWingg0

using Plots, ColorSchemes

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
labels = ["\$k_2 = $(k2)\$" for k2 in k2Range]
gr()

# Root locus
if λ == 1
    dampLim = [-10,2]
    freqLim = [0,50]
elseif λ == 2
    dampLim = [-10,2]
    freqLim = [0,80]
end
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        plot!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
if λ == 1
    OOPtext = "OOP"
    TIP1text = "1st T-IP"
    TIP2text = "2nd T-IP"
    Ttext = "T"
    IPtext = "IP"
    OOPtextPos = [-7.5, 20]
    TIP1textPos = [1, 12]
    TIP2textPos = [0.25, 33]
    annotate!(OOPtextPos[1], OOPtextPos[2], text(OOPtext, tsz))
    annotate!(TIP1textPos[1], TIP1textPos[2], text(TIP1text, tsz))
    annotate!(TIP2textPos[1], TIP2textPos[2], text(TIP2text, tsz))
    quiver!([OOPtextPos[1],OOPtextPos[1]+0.5,OOPtextPos[1]+0.5], [OOPtextPos[2]-2,OOPtextPos[2]-0.5,OOPtextPos[2]+1], quiver=([0,1.75,2.0], [-10,-3.5,8]), arrow=:closed, linecolor=:black)
end
display(plt_RL)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_rootlocus_lambda",λ,".pdf"))

# Root locus - focus on 1st T-IP mode
if λ == 1
    dampLim = [-1.5,1.5]
    freqLim = [0,22]
elseif λ == 2
    dampLim = [-2,1.5]
    freqLim = [0,50]
end
plt_RL_TIP1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    for mode in 1:nModes
        plot!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL_TIP1)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_rootlocus_TIP1_lambda",λ,".pdf"))

# Root locus - focus on 2nd T-IP mode
if λ == 1
    dampLim = [-1.5,2]
    freqLim = [20,50]
elseif λ == 2
    dampLim = [-3,2]
    freqLim = [20,80]
end
plt_RL_TIP2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    for mode in 1:nModes
        plot!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL_TIP2)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_rootlocus_TIP2_lambda",λ,".pdf"))

# V-g-f
if λ == 1
    dampLim = [-0.2,0.15]
    freqLim = [0,50]
elseif λ == 2
    dampLim = [-0.2,0.15]
    freqLim = [0,80]
end
for (i,k2) in enumerate(k2Range)
    plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=freqLim, tickfont=font(ts), guidefont=font(12))
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=dampLim, yticks=-0.5:.05:0.5, tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
    plot!(URange,zeros(length(URange)), c=:gray, lw=lw, ls=:dash, label=false)
    for mode in 1:nModes
        plot!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    if λ == 1
        if k2 == -0.015
            annotate!(plt_Vg, 34, 0.05, text("1st T-IP", tsz))
            annotate!(plt_Vg, 45, 0.075, text("2nd T-IP", tsz))
            quiver!([34], [0.03], quiver=([2.8], [-0.06]), arrow=:closed, linecolor=:black)
        elseif k2 == 0.0
            annotate!(plt_Vg, 34, 0.05, text("1st T-IP", tsz))
            annotate!(plt_Vg, 46, 0.075, text("2nd T-IP", tsz))
            quiver!([34], [0.03], quiver=([3], [-0.07]), arrow=:closed, linecolor=:black)
        elseif k2 == 0.015
            annotate!(plt_Vg, 41, 0.03, text("1st T-IP", tsz))
            annotate!(plt_Vg, 47, 0.075, text("2nd T-IP", tsz))
        elseif k2 == 0.03
            annotate!(plt_Vg, 42.5, 0.05, text("1st T-IP", tsz))
            annotate!(plt_Vg, 47, 0.075, text("2nd T-IP", tsz))
        elseif k2 == 0.045
            annotate!(plt_Vg, 22, 0.1, text("1st T-IP", tsz))
            annotate!(plt_Vg, 47, 0.075, text("2nd T-IP", tsz))
            quiver!([22], [0.1-0.02], quiver=([0], [-0.06]), arrow=:closed, linecolor=:black)
        end
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_Vgf",i,"_lambda",λ,".pdf"))
end

# V-g-f: colored by mode
if λ == 1
    modeColors = [:blue, :orange, :green, :red, :purple]
    modeLabels = ["OOP1" "T-IP1" "OOP2" "OOP3" "T-IP2"]
    modeColorOrder = [vcat(1:nModes), vcat(1:nModes), vcat(1:nModes), vcat(1,2,3,5,4), vcat(1,3,2,5,4)]
    for (i,k2) in enumerate(k2Range)
        plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50], tickfont=font(ts), guidefont=font(12))
        for mode in 1:nModes
            plot!(URange, modeFrequencies[i,mode], c=modeColors[modeColorOrder[i][mode]], lw=lw, label=false)
        end
        plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.2,0.2], yticks=-0.2:.05:0.2, tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=(0.15,1.1))
        plot!(URange,zeros(length(URange)), c=:gray, lw=lw, ls=:dash, label=false)
        for mode in 1:nModes
            if i==2
                plot!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=modeColors[modeColorOrder[i][mode]], lw=lw, label=modeLabels[mode])
            else
                plot!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=modeColors[modeColorOrder[i][mode]], lw=lw, label=false)
            end
        end
        plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
        display(plt_Vgf)
        savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_Vgf",i,"_byMode_lambda",λ,".pdf"))
    end
end

# Flutter onset speeds vs k2
if λ == 1
    ULim = [0,45]
elseif λ == 2
    ULim = [0,65]
end
plt_Uf = plot(xlabel="\$k_2\$ [1/m]", ylabel="Flutter speed [m/s]", ylims=ULim, xticks=k2Range, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomleft)
plot!(k2Range,flutterOnsetSpeedWingg0, marker=(:circle, 10), msw=msw, c=:black, lw=lw, ls=:dash, label="Wing - zero pitch, no gravity")
plot!(k2Range,flutterOnsetSpeedWingMatched, marker=(:square, 10), msw=msw, c=:black, lw=lw, ls=:solid, label="Wing - trim-matched pitch")
display(plt_Uf)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_speedOn_lambda",λ,".pdf"))

# Flutter onset speeds vs tip OOP position
plt_Uf = plot(xlabel="Tip OOP position [% semispan]", ylabel="Flutter speed [m/s]", xlims=[-40,100], ylims=[0,80], xticks=vcat(-40:10:100), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
plot!(flutterOnsetTipOOPWingg0,flutterOnsetSpeedWingg0, marker=(:circle, 5), msw=msw, c=:black, lw=lw, ls=:dash, label="Wing - zero pitch, no gravity")
plot!(flutterOnsetTipOOPWingMatched,flutterOnsetSpeedWingMatched, marker=(:square, 5), msw=msw, c=:black, lw=lw, ls=:solid, label="Wing - trim matched pitch")
display(plt_Uf)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_speedOn_OOP_lambda",λ,".pdf"))

println("Finished cHALEwing_flutter_matchedPitch_k2_range.jl")