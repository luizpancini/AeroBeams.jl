using AeroBeams

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Altitude
h = 20e3

# Discretization
nElemWing = 40
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Set number of vibration modes
nModes = 5

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
URange = collect(20:0.5:45)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))

trimAoA = fill(NaN, length(k2Range), length(URange))

untrackedFreqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
damps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeDampings = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeFrequencies = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
flutterOnsetSpeedOfMode = fill(NaN, length(k2Range), nModes)
flutterOffsetSpeedOfMode = fill(NaN, length(k2Range), nModes)
flutterOnsetSpeed = fill(NaN, length(k2Range))

highestConvUindex = Array{Int64}(undef,length(k2Range))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        # Create and trim problem
        trimProblem[i,j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        println("Trim AoA = $(trimAoA[i,j]*180/π)")
        # Model for eigen problem
        wingModel,_ = create_SMW(aeroSolver=aeroSolver,airspeed=U,nElem=nElemWing,altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA[i,j])
        # Set initial guess solution as previous known solution
        x0Eig = j == 1 ? zeros(0) : eigenProblem[i,j-1].x
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[1,Inf64],systemSolver=NReigen,x0=x0Eig)
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
    end
    # Frequencies and dampings after mode tracking
    if modeTracking
        freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
    end
    # Flutter speed of each mode
    for mode in 1:nModes
        iOnset = findfirst(j -> modeDampings[i,mode][j] < 0 && modeDampings[i,mode][j+1] > 0, 1:length(URange)-1)
        iOffset = findfirst(j -> modeDampings[i,mode][j] > 0 && modeDampings[i,mode][j+1] < 0, 1:length(URange)-1)
        if modeDampings[i,mode][1] > 0
            flutterOnsetSpeedOfMode[i,mode] = URange[1]
        elseif isnothing(iOnset)
            flutterOnsetSpeedOfMode[i,mode] = Inf64
        else
            flutterOnsetSpeedOfMode[i,mode] = interpolate(modeDampings[i,mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
        end
        if isnothing(iOffset) || isnothing(iOnset)
            flutterOffsetSpeedOfMode[i,mode] = Inf64
        else
            flutterOffsetSpeedOfMode[i,mode] = interpolate(-modeDampings[i,mode][iOffset:iOffset+1],URange[iOffset:iOffset+1],0)
        end
    end
    flutterOnsetSpeed[i] = minimum(filter(!isinf,flutterOnsetSpeedOfMode[i,:]),init=Inf64)
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/outputs/figures/cHALEwing_flutter_matchedPitch_k2_range.jl"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], range(0, 1, length(k2Range)))
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
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,1.5], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
OOPtext = "OOP"
TIP1text = "1st T-IP"
TIP2text = "2nd T-IP"
Ttext = "T"
IPtext = "IP"
OOPtextPos = [-7.5, 20]
TIP1textPos = [1, 10]
TIP2textPos = [0.25, 33]
annotate!(OOPtextPos[1], OOPtextPos[2], text(OOPtext, tsz))
annotate!(TIP1textPos[1], TIP1textPos[2], text(TIP1text, tsz))
annotate!(TIP2textPos[1], TIP2textPos[2], text(TIP2text, tsz))
quiver!([OOPtextPos[1],OOPtextPos[1]+0.5,OOPtextPos[1]+0.5], [OOPtextPos[2]-2,OOPtextPos[2]-0.5,OOPtextPos[2]+1], quiver=([0,1.75,3.0], [-10,-3.5,13]), arrow=:closed, linecolor=:black)
display(plt_RL)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_rootlocus.pdf"))

# V-g-f
Vgf_k2_ind = [2,5]
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,41], tickfont=font(ts), guidefont=font(12))
for (i,k2) in enumerate(k2Range)
    if !(i in Vgf_k2_ind)
        continue
    end
    for mode in 1:nModes
        scatter!(URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.2,0.1], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for (i,k2) in enumerate(k2Range)
    if !(i in Vgf_k2_ind)
        continue
    end
    for mode in 1:nModes
        scatter!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_Vgf.pdf"))

# Flutter onset speeds vs k2
plt_Uf = plot(xlabel="\$k_2\$", ylabel="Flutter speed [m/s]", ylims=[0,45], xticks=k2Range, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomleft)
plot!(k2Range,[22.42; 32.74; 22.42; 16.97; 12.50], marker=(:circle, 10), msw=0, c=:black, lw=lw, ls=:dash, label="Wing - zero pitch, no gravity")
plot!(k2Range,flutterOnsetSpeed, marker=(:square, 10), msw=0, c=:black, lw=lw, ls=:solid, label="Wing - trim matched pitch")
display(plt_Uf)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_speedOn.pdf"))

println("Finished cHALEwing_flutter_matchedPitch_k2_range.jl")