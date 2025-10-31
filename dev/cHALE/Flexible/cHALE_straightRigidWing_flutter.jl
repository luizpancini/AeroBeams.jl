using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1e2

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
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,pseudoInverseMethod=:dampedLeastSquares,displayStatus=false)

# Set number of vibration modes
nModes = 6

# Airspeed range
URange = vcat(20:1:35)

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
μu = 1e-1
μp = 1e-1
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

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_straightRigidWing_flutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at U = 25 m/s
iU25 = findfirst(x->x≈25,URange)
plt_modes = plot_mode_shapes(eigenProblem[iU25],nModes=6,scale=5,view=(30,30),modalColorScheme=:rainbow,legendPos=:top,plotLimits=([-20,-20,-20],[20,20,20]),ΔuDef=[0,0,-5],save=true,savePath=string(relPath,string("/cHALE_flutter_k2_",k2,"_modeShapes.pdf")))
display(plt_modes)

# Plot configurations
modeColors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
gr()

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,0.1], ylims=[0,6], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
display(plt_RL)
savefig(string(absPath,"/cHALE_straightRigidWing_flutter_k2_",k2,"_rootlocus.pdf"))

# Root locus - phugoid focus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.06,0], ylims=[0,0.5], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
display(plt_RL)
savefig(string(absPath,"/cHALE_straightRigidWing_flutter_k2_",k2,"_rootlocus_phugoid.pdf"))

# V-g-f
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,6], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-1.5,0.05], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALE_straightRigidWing_flutter_k2_",k2,"_Vgf.pdf"))

println("Finished cHALE_straightRigidWing_flutter.jl")