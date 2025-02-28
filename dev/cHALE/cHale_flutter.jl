using AeroBeams

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Bending pre-curvature
k2 = 0.030

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
maxIter = 30
σ0 = 1
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Set number of vibration modes
nModes = 20

# Airspeed range
URange = vcat(20:1:45)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(URange))
eigenProblem = Array{EigenProblem}(undef,length(URange))

untrackedFreqs = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for U in 1:length(URange)]
damps = [fill(NaN64, nModes) for U in 1:length(URange)]
modeDampings = [fill(NaN64, nModes) for U in 1:length(URange)]
modeFrequencies = [fill(NaN64, nModes) for U in 1:length(URange)]

# Set attachment springs
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=1e-2*[1; 1; 1],kp=1e-2*[1; 1; 1])
spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=1e-2*[1; 1; 1],kp=1e-2*[1; 1; 1])

# Sweep airspeed
for (j,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
    # Add springs
    add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
    # Update model
    cHALEtrim.skipValidationMotionBasisA = true
    update_model!(cHALEtrim)
    # Set initial guess solution as previous known solution
    x0Trim = j == 1 ? zeros(0) : trimProblem[j-1].x
    # Create and trim problem
    trimProblem[j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
    solve!(trimProblem[j])
    # Skip if unconverged
    if !trimProblem[j].systemSolver.convergedFinalSolution
        break
    end
    # Extract trim variables
    trimAoA = trimProblem[j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
    trimThrust = stabilizersAero ? trimProblem[j].x[end-1]*trimProblem[j].model.forceScaling : trimProblem[j].x[end]*trimProblem[j].model.forceScaling
    trimδ = stabilizersAero ? trimProblem[j].x[end] : 0
    println("Trim AoA = $(trimAoA*180/π), trim thrust = $(trimThrust), trim δ = $(trimδ*180/π)")
    # Model for eigen problem
    cHALEeigen,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ,thrust=trimThrust,k2=k2,hasInducedDrag=hasInducedDrag)
    # Create and solve eigen problem
    eigenProblem[j] = create_EigenProblem(model=cHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],jacobian=trimProblem[j].jacobian[1:end,1:end-trimProblem[j].model.nTrimVariables],inertia=trimProblem[j].inertia,refTrimProblem=trimProblem[j])
    solve_eigen!(eigenProblem[j])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[j] = eigenProblem[j].frequenciesOscillatory
    untrackedDamps[j] = round_off!(eigenProblem[j].dampingsOscillatory,1e-8)
    untrackedEigenvectors[j] = eigenProblem[j].eigenvectorsOscillatoryCplx
end

# Frequencies and dampings after mode tracking
if modeTracking
    freqs,damps,_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and dampings by mode
modeDampings = [fill(NaN64, nModes) for U in 1:length(URange)]
modeFrequencies = [fill(NaN64, nModes) for U in 1:length(URange)]
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[j][mode] for j in eachindex(URange)]
    modeDampings[mode] = [damps[j][mode] for j in eachindex(URange)]
end

# Flutter speed of each mode
flutterOnsetSpeedOfMode = fill(NaN, nModes)
flutterOffsetSpeedOfMode = fill(NaN, nModes)
for mode in 1:nModes
    iOnset = findfirst(j -> modeDampings[mode][j] < 0 && modeDampings[mode][j+1] > 0, 1:length(URange)-1)
    iOffset = findfirst(j -> modeDampings[mode][j] > 0 && modeDampings[mode][j+1] < 0, 1:length(URange)-1)
    if isnothing(iOnset)
        flutterOnsetSpeedOfMode[mode] = Inf64
    else
        flutterOnsetSpeedOfMode[mode] = interpolate(modeDampings[mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
    end
    if isnothing(iOffset) || isnothing(iOnset)
        flutterOffsetSpeedOfMode[mode] = Inf64
    else
        flutterOffsetSpeedOfMode[mode] = interpolate(-modeDampings[mode][iOffset:iOffset+1],URange[iOffset:iOffset+1],0)
    end
end
flutterOnsetSpeed = minimum(filter(!isinf,flutterOnsetSpeedOfMode),init=Inf64)

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/outputs/figures/cHale_flutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at lowest airspeed
plt_modes = plot_mode_shapes(eigenProblem[1],nModes=6,scale=5,view=(30,30),modalColorScheme=:rainbow,modeLabels=["Phugoid","1st LD","1st Sym OOP","2nd LD","2nd Sym OOP","1st Asym OOP"],legendPos=:top,plotLimits=([-20,-20,-20],[20,20,20]),ΔuDef=[-5,10,25],save=true,savePath=string(relPath,string("/cHale_flutter_k2_",k2,"_modeShapes.pdf")))
display(plt_modes)

# Plot configurations
colors = cgrad(:thermal, length(URange), categorical=true)
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
gr()

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-4,1.5], ylims=[5,21], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
display(plt_RL)
savefig(string(absPath,"/cHale_flutter_k2_",k2,"_rootlocus.pdf"))

println("Finished cHale_flutter.jl")