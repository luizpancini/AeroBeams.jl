using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Bending pre-curvature
k2 = 0.0

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
    nElemWing = 40
elseif λ > 1
    nElemWing = 40
end
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,pseudoInverseMethod=:dampedLeastSquares,displayStatus=true)

# Set number of vibration modes
nModes = 10

# Airspeed range
if λ == 1
    URange = [20]
elseif λ == 2
    URange = [20]
else
    URange = [20]
end

# Damping ratio tolerance for flutter detection
σtol = 0e-2

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(URange))
eigenProblem = Array{EigenProblem}(undef,length(URange))

untrackedFreqs = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for U in 1:length(URange)]
damps = [fill(NaN64, nModes) for U in 1:length(URange)]
modeDampings = [fill(NaN64, length(URange)) for mode in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for mode in 1:nModes]

# Set attachment springs
μu = 1e-2
μp = 0e-2
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])
spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    cHALEtrim,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
    # Set initial guess solution as previous known solution
    x0Trim = i == 1 ? zeros(0) : trimProblem[i-1].x
    # Create and trim problem
    trimProblem[i] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
    solve!(trimProblem[i])
    # Skip if unconverged
    if !trimProblem[i].systemSolver.convergedFinalSolution
        break
    end
    # Extract trim variables
    trimAoA = trimProblem[i].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
    trimThrust = stabilizersAero ? trimProblem[i].x[end-1]*trimProblem[i].model.forceScaling : trimProblem[i].x[end]*trimProblem[i].model.forceScaling
    trimδ = stabilizersAero ? trimProblem[i].x[end] : 0
    println("Trim AoA = $(trimAoA*180/π), trim thrust = $(trimThrust), trim δ = $(trimδ*180/π)")
    # Model for trim problem with springs
    cHALEtrimSpringed,_,_,tailBoomSpringed,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
    # Add springs
    add_springs_to_beam!(beam=tailBoomSpringed,springs=[spring1,spring2])
    # Update model
    update_model!(cHALEtrimSpringed)
    # Create and solve trim problem with springs
    trimProblemSpringed = create_TrimProblem(model=cHALEtrimSpringed,systemSolver=NRtrim,x0=trimProblem[i].x)
    solve!(trimProblemSpringed)
    # Retrieve and compare trim outputs
    trimAoASpringed = (trimProblemSpringed.aeroVariablesOverσ[end][cHALEtrimSpringed.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblemSpringed.aeroVariablesOverσ[end][cHALEtrimSpringed.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
    trimThrustSpringed = trimProblemSpringed.x[end-1]*trimProblemSpringed.model.forceScaling
    trimδSpringed = trimProblemSpringed.x[end]
    println("Trim outputs ratios springed/nominal: AoA = $(trimAoASpringed/trimAoA), T = $(trimThrustSpringed/trimThrust), δ = $(trimδSpringed/trimδ)")
    # Model for eigen problem
    cHALEeigen,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ,thrust=trimThrust,k2=k2,hasInducedDrag=hasInducedDrag)
    # Create and solve eigen problem
    eigenProblem[i] = create_EigenProblem(model=cHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-1,Inf],refTrimProblem=trimProblemSpringed)
    solve_eigen!(eigenProblem[i])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem[i].frequenciesOscillatory
    untrackedDamps[i] = eigenProblem[i].dampingsOscillatory
    untrackedEigenvectors[i] = eigenProblem[i].eigenvectorsOscillatoryCplx
end

# Frequencies and dampings after mode tracking
freqs,damps,_ = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and dampings by mode
modeDampings = [fill(NaN64, length(URange)) for m in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for m in 1:nModes]
modeDampingRatios = [fill(NaN64, length(URange)) for m in 1:nModes]
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

# Flutter onset/offset data of each mode
flutterOnsetMode = fill(0,0)
flutterOffsetMode = fill(0,0)
flutterOnsetSpeedOfMode = [fill(0.0,0) for _ in 1:nModes]
flutterOnsetFreqOfMode = [fill(0.0,0) for _ in 1:nModes]
flutterOffsetSpeedOfMode = [fill(0.0,0) for _ in 1:nModes]
flutterOffsetFreqOfMode = [fill(0.0,0) for _ in 1:nModes]
for mode in 1:nModes
    iOnset = findall(j -> modeDampings[mode][j] < σtol && modeDampings[mode][j+1] > σtol, 1:length(URange)-1)
    iOffset = findall(j -> modeDampings[mode][j] > σtol && modeDampings[mode][j+1] < σtol, 1:length(URange)-1)
    if modeDampingRatios[mode][1] > σtol
        push!(flutterOnsetMode,mode)
        push!(flutterOnsetSpeedOfMode[mode],URange[1])
        push!(flutterOnsetFreqOfMode[mode],modeFrequencies[mode][1])
    end
    if !isempty(iOnset)
        for iO in iOnset
            push!(flutterOnsetMode,mode)
            push!(flutterOnsetSpeedOfMode[mode],interpolate(modeDampings[mode][iO:iO+1],URange[iO:iO+1],σtol))
            push!(flutterOnsetFreqOfMode[mode],interpolate(modeDampings[mode][iO:iO+1],modeFrequencies[mode][iO:iO+1],σtol))
        end
    end
    if !(isempty(iOffset) || isempty(iOnset))
        for iO in iOffset
            push!(flutterOffsetMode,mode)
            push!(flutterOffsetSpeedOfMode[mode],interpolate(-modeDampings[mode][iO:iO+1],URange[iO:iO+1],-σtol))
            push!(flutterOffsetFreqOfMode[mode],interpolate(-modeDampings[mode][iO:iO+1],modeFrequencies[mode][iO:iO+1],-σtol))
        end
    end
end

# All flutter onset/offset speeds and frequencies
flutterOnsetSpeedsAll = vcat(filter(!isempty,flutterOnsetSpeedOfMode)...)
flutterOnsetFreqsAll = vcat(filter(!isempty,flutterOnsetFreqOfMode)...)
flutterOffsetSpeedsAll = vcat(filter(!isempty,flutterOffsetSpeedOfMode)...)
flutterOffsetFreqsAll = vcat(filter(!isempty,flutterOffsetFreqOfMode)...)

# Lowest flutter onset speed and corresponding frequency
if !isempty(flutterOnsetSpeedsAll)
    iLowest = sortperm(flutterOnsetSpeedsAll)[1]
    flutterOnsetSpeed = flutterOnsetSpeedsAll[iLowest]
    flutterOnsetFreq = flutterOnsetFreqsAll[iLowest]
end

using Plots, ColorSchemes

colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_flutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# # Mode shapes at lowest airspeed
# plt_modes = plot_mode_shapes(eigenProblem[1],nModes=6,scale=5,view=(30,30),modalColorScheme=colors,modeLabels=["Phugoid","1st LD","1st Sym OOP","2nd LD","2nd Sym OOP","1st Asym OOP"],legendPos=:top,plotLimits=([-20,-20,-20],[20,20,20]),ΔuDef=[-5,10,25],save=true,savePath=string(relPath,string("/cHALE_flutter_k2_",k2,"_modeShapes_lambda",λ,".pdf")))
# display(plt_modes)

# Modes animation
modesAnim = plot_mode_shapes_animation(eigenProblem[1],scale=10,view=(60,15),modalColorScheme=colors,legendPos=:top,modes2plot=[3,5],nFramesPerCycle=41,plotSteady=false,plotBCs=false,plotAxes=false,legendFontSize=10,fps=20,save=true,savePath=string(relPath,"/cHALE_flutter_modes_k2_",k2,"_lambda",λ,".gif"),displayProgress=true)
display(modesAnim)

# Plot configurations
modeColors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold], nModes)
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
gr()

# Root locus
if λ == 1
    dampLim = [-10,2]
    freqLim = [0,50]
elseif λ == 2
    dampLim = [-10,2]
    freqLim = [0,80]
else
    dampLim = [-10,2]
    freqLim = [0,100]
end
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
for mode in 1:nModes
    plot!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
    plot!([modeDampings[mode][1]], [modeFrequencies[mode][1]], c=modeColors[mode], shape=:circle, ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RL)
savefig(string(absPath,"/cHALE_flutter_lambda",λ,"_k2_",k2,"_rootlocus.pdf"))

# Root locus - low frequency modes
if λ == 1
    dampLim = [-6,2]
    freqLim = [0,5]
elseif λ == 2
    dampLim = [-6,2]
    freqLim = [0,10]
end
plt_RLlow = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
for mode in 1:nModes
    plot!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
    plot!([modeDampings[mode][1]], [modeFrequencies[mode][1]], c=modeColors[mode], shape=:circle, ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RLlow)
savefig(string(absPath,"/cHALE_flutter_lambda",λ,"_k2_",k2,"_rootlocus_low.pdf"))

# Root locus (zoom on T-IP modes)
if λ == 1
    dampLim = [-2.5,2]
    freqLim = [2,45]
elseif λ == 2
    dampLim = [-2.5,4]
    freqLim = [5,80]
end     
plt_RLTIP = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    plot!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
    plot!([modeDampings[mode][1]], [modeFrequencies[mode][1]], c=modeColors[mode], shape=:circle, ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RLTIP)
savefig(string(absPath,"/cHALE_flutter_lambda",λ,"_k2_",k2,"_rootlocus_TIP.pdf"))

# Root locus - phugoid
if λ == 1
    dampLim = [-0.15,0.15]
    freqLim = [0,0.5]
elseif λ == 2
    dampLim = [-0.15,0.35]
    freqLim = [0,0.5]
end
plt_RLphu = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=dampLim, ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
for mode in 1:nModes
    plot!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
    plot!([modeDampings[mode][1]], [modeFrequencies[mode][1]], c=modeColors[mode], shape=:circle, ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RLphu)
savefig(string(absPath,"/cHALE_flutter_lambda",λ,"_k2_",k2,"_rootlocus_phugoid.pdf"))

# V-g-f
if λ == 1
    dampLim = [-0.3,0.2]
    freqLim = [0,50]
elseif λ == 2
    dampLim = [-0.3,0.25]
    freqLim = [0,80]
end
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=freqLim, tickfont=font(ts), guidefont=font(12))
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=dampLim, tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
plot!(plt_Vg,URange,zeros(length(URange)), c=:gray, lw=lw, ls=:dash, label=false)
for mode in 1:nModes
    plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALE_flutter_lambda",λ,"_k2_",k2,"_Vgf.pdf"))

println("Finished cHALE_flutter.jl")