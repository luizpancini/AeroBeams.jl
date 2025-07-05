using AeroBeams, LinearInterpolations, JLD2

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

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
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
relTol = 1e-8
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,relativeTolerance=relTol,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)
NReigenWing = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=0.5,displayStatus=false)

# Set number of vibration modes
nModesAircraft = 30
nModesWing = 20

# Bending pre-curvature
k2 = 0.045

# Damping ratio tolerance for flutter detection at lower airspeed range
σtol = 2e-3

# Airspeed range
if λ == 1
    URange = unique(sort(vcat(20:0.5:50,42.3,42.4)))
elseif λ == 2
    URange = unique(sort(vcat(20:0.5:65,46.5:0.25:47.5,60:0.25:65)))
end

# Initialize aircraft outputs
trimProblem = Array{TrimProblem}(undef,length(URange))
trimAoA = fill(NaN, length(URange))
trimThrust = fill(NaN, length(URange))
trimδ = fill(NaN, length(URange))
eigenProblemAircraft = Array{EigenProblem}(undef,length(URange))
untrackedFreqsAircraft = [fill(NaN64, nModesAircraft) for U in 1:length(URange)]
untrackedDampsAircraft = [fill(NaN64, nModesAircraft) for U in 1:length(URange)]
untrackedEigenvectorsAircraft = [fill(NaN64+im*NaN64, nModesAircraft, nModesAircraft) for U in 1:length(URange)]
freqsAircraft = [fill(NaN64, nModesAircraft) for U in 1:length(URange)]
dampsAircraft = [fill(NaN64, nModesAircraft) for U in 1:length(URange)]
modeDampingsAircraft = [fill(NaN64, length(URange)) for mode in 1:nModesAircraft]
modeDampingRatiosAircraft = [fill(NaN64, length(URange)) for mode in 1:nModesAircraft]
modeFrequenciesAircraft = [fill(NaN64, length(URange)) for mode in 1:nModesAircraft]
flutterSpeedOfModeAircraft = fill(NaN, nModesAircraft)
flutterFreqOfModeAircraft = fill(NaN, nModesAircraft)
flutterTipOOPOfModeAircraft = fill(NaN, nModesAircraft)

# Initialize trimmed wing outputs
eigenProblemWing = Array{EigenProblem}(undef,length(URange))
untrackedFreqsWing = [fill(NaN64, nModesWing) for U in 1:length(URange)]
untrackedDampsWing = [fill(NaN64, nModesWing) for U in 1:length(URange)]
untrackedEigenvectorsWing = [fill(NaN64+im*NaN64, nModesWing, nModesWing) for U in 1:length(URange)]
freqsWing = [fill(NaN64, nModesWing) for U in 1:length(URange)]
dampsWing = [fill(NaN64, nModesWing) for U in 1:length(URange)]
modeDampingsWing = [fill(NaN64, length(URange)) for mode in 1:nModesWing]
modeDampingRatiosWing = [fill(NaN64, length(URange)) for mode in 1:nModesWing]
modeFrequenciesWing = [fill(NaN64, length(URange)) for mode in 1:nModesWing]
flutterSpeedOfModeWing = fill(NaN, nModesWing)
flutterFreqOfModeWing = fill(NaN, nModesWing)
flutterTipOOPOfModeWing = fill(NaN, nModesWing)

# Initialize common outputs
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))

# ELement ranges
elemRangeRightWing = 1 + div(nElemWing,2) : nElemWing

# Set attachment springs
μu = 1e-2
μp = 1e-2
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])
spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])

# Undeformed jig-shape properties
model,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=URange[1],nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,k2=k2,hasInducedDrag=hasInducedDrag)
x1_0 = vcat([vcat(model.elements[e].r_n1[1],model.elements[e].r_n2[1]) for e in elemRangeRightWing]...)
x3_0 = vcat([vcat(model.elements[e].r_n1[3],model.elements[e].r_n2[3]) for e in elemRangeRightWing]...)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
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
    trimAoA[i] = trimProblem[i].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
    trimThrust[i] = stabilizersAero ? trimProblem[i].x[end-1]*trimProblem[i].model.forceScaling : trimProblem[i].x[end]*trimProblem[i].model.forceScaling
    trimδ[i] = stabilizersAero ? trimProblem[i].x[end] : 0
    println("Trim AoA = $(trimAoA[i]*180/π), trim thrust = $(trimThrust[i]), trim δ = $(trimδ[i]*180/π)")
    # Displacements over span
    u1_of_x1 = vcat([vcat(trimProblem[i].nodalStatesOverσ[end][e].u_n1[1],trimProblem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in elemRangeRightWing]...)
    u3_of_x1 = vcat([vcat(trimProblem[i].nodalStatesOverσ[end][e].u_n1[3],trimProblem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in elemRangeRightWing]...)
    u1_of_x1 .-= u1_of_x1[1]
    u3_of_x1 .-= u3_of_x1[1]
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1_of_x1
    x3_def[i] = x3_0 .+ u3_of_x1
    # Model for trim problem with springs
    cHALEtrimSpringed,_,_,tailBoomSpringed,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
    # Add springs
    add_springs_to_beam!(beam=tailBoomSpringed,springs=[spring1,spring2])
    # Update model
    cHALEtrimSpringed.skipValidationMotionBasisA = true
    update_model!(cHALEtrimSpringed)
    # Create and solve trim problem with springs
    trimProblemSpringed = create_TrimProblem(model=cHALEtrimSpringed,systemSolver=NRtrim,x0=trimProblem[i].x)
    solve!(trimProblemSpringed)
    # Retrieve and compare trim outputs
    trimAoASpringed = (trimProblemSpringed.aeroVariablesOverσ[end][cHALEtrimSpringed.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblemSpringed.aeroVariablesOverσ[end][cHALEtrimSpringed.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
    trimThrustSpringed = trimProblemSpringed.x[end-1]*trimProblemSpringed.model.forceScaling
    trimδSpringed = trimProblemSpringed.x[end]
    println("Trim outputs ratios springed/nominal: AoA = $(trimAoASpringed/trimAoA[i]), T = $(trimThrustSpringed/trimThrust[i]), δ = $(trimδSpringed/trimδ[i])")
    # Model for eigen problem
    cHALEeigen,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδSpringed,thrust=trimThrustSpringed,k2=k2,hasInducedDrag=hasInducedDrag)
    # Create and solve eigen problem
    eigenProblemAircraft[i] = create_EigenProblem(model=cHALEeigen,nModes=nModesAircraft,frequencyFilterLimits=[1e-1,Inf],refTrimProblem=trimProblemSpringed)
    solve_eigen!(eigenProblemAircraft[i])
    # Frequencies, dampings and eigenvectors of aircraft
    untrackedFreqsAircraft[i] = eigenProblemAircraft[i].frequenciesOscillatory
    untrackedDampsAircraft[i] = round_off!(eigenProblemAircraft[i].dampingsOscillatory,1e-8)
    untrackedEigenvectorsAircraft[i] = eigenProblemAircraft[i].eigenvectorsOscillatoryCplx
    # Model for trimmed wing eigen problem
    wingModel,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA[i])
    # Set initial guess solution as previous known solution
    x0Eig = i == 1 ? zeros(0) : eigenProblemWing[i-1].x
    # Create and solve eigen problem for trimmed wing
    eigenProblemWing[i] = create_EigenProblem(model=wingModel,nModes=nModesWing,frequencyFilterLimits=[1e-2,Inf],systemSolver=NReigenWing,x0=x0Eig)
    solve!(eigenProblemWing[i])
    # Frequencies, dampings and eigenvectors of trimmed wing
    untrackedFreqsWing[i] = eigenProblemWing[i].frequenciesOscillatory
    untrackedDampsWing[i] = round_off!(eigenProblemWing[i].dampingsOscillatory,1e-8)
    untrackedEigenvectorsWing[i] = eigenProblemWing[i].eigenvectorsOscillatoryCplx
end

# Aircraft frequencies and dampings after mode tracking
freqsAircraft,dampsAircraft,_ = mode_tracking_hungarian(URange,untrackedFreqsAircraft,untrackedDampsAircraft,untrackedEigenvectorsAircraft)

# Trimmed wing frequencies and dampings after mode tracking
freqsWing,dampsWing,_ = mode_tracking_hungarian(URange,untrackedFreqsWing,untrackedDampsWing,untrackedEigenvectorsWing)

# Separate frequencies and dampings by mode
for mode in 1:nModesAircraft
    modeFrequenciesAircraft[mode] = [freqsAircraft[i][mode] for i in eachindex(URange)]
    modeDampingsAircraft[mode] = [dampsAircraft[i][mode] for i in eachindex(URange)]
    modeDampingRatiosAircraft[mode] = modeDampingsAircraft[mode]./modeFrequenciesAircraft[mode]
end
for mode in 1:nModesWing
    modeFrequenciesWing[mode] = [freqsWing[i][mode] for i in eachindex(URange)]
    modeDampingsWing[mode] = [dampsWing[i][mode] for i in eachindex(URange)]
    modeDampingRatiosWing[mode] = modeDampingsWing[mode]./modeFrequenciesWing[mode]
end

# Tolerance for true damping ratio crossing
ϵ = -5e-4

# Flutter speed of each mode for aircraft
for mode in 1:nModesAircraft
    iOnset = findfirst(i -> modeDampingsAircraft[mode][i] < 0 && modeDampingsAircraft[mode][i+1] > 0, 1:length(URange)-1)
    if modeDampingsAircraft[mode][1]/modeFrequenciesAircraft[mode][1] > σtol
        flutterSpeedOfModeAircraft[mode] = URange[1]
        flutterFreqOfModeAircraft[mode] = modeFrequenciesAircraft[mode][1]
        flutterTipOOPOfModeAircraft[mode] = x3_def[1][end]/16*100
        println("Aircraft, mode $mode flutter: i=$iOnset, U=$(flutterSpeedOfModeAircraft[mode]), ω=$(flutterFreqOfModeAircraft[mode])")
    elseif isnothing(iOnset)
        flutterSpeedOfModeAircraft[mode] = Inf
        flutterFreqOfModeAircraft[mode] = NaN
        flutterTipOOPOfModeAircraft[mode] = NaN
    else
        # Check false crossing (mode swapping)
        if modeDampingRatiosAircraft[mode][iOnset] * modeDampingRatiosAircraft[mode][iOnset+1] < ϵ
            flutterSpeedOfModeAircraft[mode] = Inf
            flutterFreqOfModeAircraft[mode] = NaN
            flutterTipOOPOfModeAircraft[mode] = NaN
            continue
        end
        flutterSpeedOfModeAircraft[mode] = interpolate(modeDampingsAircraft[mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
        flutterFreqOfModeAircraft[mode] = interpolate(modeDampingsAircraft[mode][iOnset:iOnset+1],modeFrequenciesAircraft[mode][iOnset:iOnset+1],0)
        flutterTipOOPOfModeAircraft[mode] = (x3_def[iOnset][end]+x3_def[iOnset+1][end])/2/16*100
        println("Aircraft, mode $mode flutter: i=$iOnset, U=$(flutterSpeedOfModeAircraft[mode]), ω=$(flutterFreqOfModeAircraft[mode])")
    end
end

# All aircraft flutter onset speeds and frequencies
flutterSpeedAircraftAll = filter(!isinf,flutterSpeedOfModeAircraft)
flutterFreqAircraftAll = filter(!isinf,flutterFreqOfModeAircraft)
fluttertipOOPAircraftAll = filter(!isinf,flutterTipOOPOfModeAircraft)
indicesFlutterOnsetAircraft = []
for speed in flutterSpeedAircraftAll
    push!(indicesFlutterOnsetAircraft, -1 + findfirst(x-> x > speed, URange))
end
sort!(indicesFlutterOnsetAircraft)
# 1st flutter onset: speed index, mode, speed, frequency, and OOP position
flutterSpeedIndexAircraft = indicesFlutterOnsetAircraft[1]
flutterModeAircraft = sortperm(flutterSpeedOfModeAircraft)[1]
flutterSpeedAircraft = flutterSpeedOfModeAircraft[flutterModeAircraft]
flutterFreqAircraft = flutterFreqOfModeAircraft[flutterModeAircraft]
flutterTipOOPAircraft = flutterTipOOPOfModeAircraft[flutterModeAircraft]

# Flutter speed of each mode for trimmed wing
for mode in 1:nModesWing
    iOnset = findfirst(i -> modeDampingsWing[mode][i] < 0 && modeDampingsWing[mode][i+1] > 0, 1:length(URange)-1)
    if modeDampingsWing[mode][1]/modeFrequenciesWing[mode][1] > σtol
        flutterSpeedOfModeWing[mode] = URange[1]
        flutterFreqOfModeWing[mode] = modeFrequenciesWing[mode][1]
        flutterTipOOPOfModeWing[mode] = x3_def[1][end]/16*100
        println("Trimmed wing, mode $mode flutter: i=$iOnset, U=$(flutterSpeedOfModeWing[mode]), ω=$(flutterFreqOfModeWing[mode])")
    elseif isnothing(iOnset)
        flutterSpeedOfModeWing[mode] = Inf
        flutterFreqOfModeWing[mode] = NaN
        flutterTipOOPOfModeWing[mode] = NaN
    else
        # Check false crossing (mode swapping)
        if modeDampingRatiosWing[mode][iOnset] * modeDampingRatiosWing[mode][iOnset+1] < ϵ
            flutterSpeedOfModeWing[mode] = Inf
            flutterFreqOfModeWing[mode] = NaN
            flutterTipOOPOfModeWing[mode] = NaN
            continue
        end
        flutterSpeedOfModeWing[mode] = interpolate(modeDampingsWing[mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
        flutterFreqOfModeWing[mode] = interpolate(modeDampingsWing[mode][iOnset:iOnset+1],modeFrequenciesWing[mode][iOnset:iOnset+1],0)
        flutterTipOOPOfModeWing[mode] = (x3_def[iOnset][end]+x3_def[iOnset+1][end])/2/16*100
        println("Trimmed wing, mode $mode flutter: i=$iOnset, U=$(flutterSpeedOfModeWing[mode]), ω=$(flutterFreqOfModeWing[mode])")
    end
end

# All trimmed wing flutter onset speeds and frequencies, in order
flutterSpeedWingAll = filter(!isinf,flutterSpeedOfModeWing)
flutterFreqWingAll = filter(!isinf,flutterFreqOfModeWing)
fluttertipOOPWingAll = filter(!isinf,flutterTipOOPOfModeWing)
indicesFlutterOnsetWing = []
for speed in flutterSpeedWingAll
    push!(indicesFlutterOnsetWing, -1 + findfirst(x-> x > speed, URange))
end
sort!(indicesFlutterOnsetWing)
# 1st flutter onset: speed index, mode, speed, frequency, and OOP position
flutterSpeedIndexWing = indicesFlutterOnsetWing[1]
flutterModeWing = sortperm(flutterSpeedOfModeWing)[1]
flutterSpeedWing = flutterSpeedOfModeWing[flutterModeWing]
flutterFreqWing = flutterFreqOfModeWing[flutterModeWing]
flutterTipOOPWing = flutterTipOOPOfModeWing[flutterModeWing]

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_compare_stability"
absPath = string(pwd(),relPath)
mkpath(absPath)
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALE_compare_stability/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)

# Save data
@save absPathData*string("x1_0_lambda",λ,"_k2",k2)*".jld2" x1_0
@save absPathData*string("x3_0_lambda",λ,"_k2",k2)*".jld2" x3_0
@save absPathData*string("URange_lambda",λ,"_k2",k2)*".jld2" URange
@save absPathData*string("x1_0_lambda",λ,"_k2",k2)*".jld2" x1_0
@save absPathData*string("x3_0_lambda",λ,"_k2",k2)*".jld2" x3_0
@save absPathData*string("x3_def_lambda",λ,"_k2",k2)*".jld2" x3_def
@save absPathData*string("x1_def_lambda",λ,"_k2",k2)*".jld2" x1_def
@save absPathData*string("x3_def_lambda",λ,"_k2",k2)*".jld2" x3_def
@save absPathData*string("freqsAircraft_lambda",λ,"_k2",k2)*".jld2" freqsAircraft
@save absPathData*string("dampsAircraft_lambda",λ,"_k2",k2)*".jld2" dampsAircraft
@save absPathData*string("flutterSpeedAircraftAll_lambda",λ,"_k2",k2)*".jld2" flutterSpeedAircraftAll
@save absPathData*string("flutterFreqAircraftAll_lambda",λ,"_k2",k2)*".jld2" flutterFreqAircraftAll
@save absPathData*string("fluttertipOOPAircraftAll_lambda",λ,"_k2",k2)*".jld2" fluttertipOOPAircraftAll
@save absPathData*string("flutterSpeedIndexAircraft_lambda",λ,"_k2",k2)*".jld2" flutterSpeedIndexAircraft
@save absPathData*string("flutterModeAircraft_lambda",λ,"_k2",k2)*".jld2" flutterModeAircraft
@save absPathData*string("flutterSpeedAircraft_lambda",λ,"_k2",k2)*".jld2" flutterSpeedAircraft
@save absPathData*string("flutterFreqAircraft_lambda",λ,"_k2",k2)*".jld2" flutterFreqAircraft
@save absPathData*string("flutterTipOOPAircraft_lambda",λ,"_k2",k2)*".jld2" flutterTipOOPAircraft
@save absPathData*string("freqsWing_lambda",λ,"_k2",k2)*".jld2" freqsWing
@save absPathData*string("dampsWing_lambda",λ,"_k2",k2)*".jld2" dampsWing
@save absPathData*string("flutterSpeedWingAll_lambda",λ,"_k2",k2)*".jld2" flutterSpeedWingAll
@save absPathData*string("flutterFreqWingAll_lambda",λ,"_k2",k2)*".jld2" flutterFreqWingAll
@save absPathData*string("fluttertipOOPWingAll_lambda",λ,"_k2",k2)*".jld2" fluttertipOOPWingAll
@save absPathData*string("flutterSpeedIndexWing_lambda",λ,"_k2",k2)*".jld2" flutterSpeedIndexWing
@save absPathData*string("flutterModeWing_lambda",λ,"_k2",k2)*".jld2" flutterModeWing
@save absPathData*string("flutterSpeedWing_lambda",λ,"_k2",k2)*".jld2" flutterSpeedWing
@save absPathData*string("flutterFreqWing_lambda",λ,"_k2",k2)*".jld2" flutterFreqWing
@save absPathData*string("flutterTipOOPWing_lambda",λ,"_k2",k2)*".jld2" flutterTipOOPWing

# Plot configurations
using Plots, ColorSchemes
modecolorsA = cgrad(:rainbow, nModesAircraft, categorical=true)
modecolorsW = cgrad(:rainbow, nModesWing, categorical=true)
colors = cgrad(:rainbow, categorical=true)
ts = 10
fs = 16
lfs = 10
tsz = 10
lw = 2
ms = 3
msw = 0
mshape = [:circle, :square]
L = 16
gr()

# Plot aircraft V-g-fs separated by mode
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModesAircraft
    plot!(URange, modeFrequenciesAircraft[mode], c=modecolorsA[mode], shape=mshape[1], ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModesAircraft
    plot!(URange, modeDampingRatiosAircraft[mode], c=modecolorsA[mode], shape=mshape[1], ms=ms, msw=msw, label=false)
end
plt_VgfmodesAircraft = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_VgfmodesAircraft)
savefig(string(absPath,"/cHALE_compare_stability_VgfAircraft_lambda_",λ,"_k2_",k2,".pdf"))

# Plot trimmed wing V-g-fs separated by mode
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModesWing
    plot!(URange, modeFrequenciesWing[mode], c=modecolorsW[mode], shape=mshape[1], ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModesWing
    plot!(URange, modeDampingRatiosWing[mode], c=modecolorsW[mode], shape=mshape[1], ms=ms, msw=msw, label=false)
end
plt_VgfmodesWing = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_VgfmodesWing)
savefig(string(absPath,"/cHALE_compare_stability_VgfWing_lambda_",λ,"_k2_",k2,".pdf"))

# Root locus - damping ratio vs frequency
plt_RLdampr = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-0.5,0.2], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label="Aircraft")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label="Trimmed wing")
for mode in 1:nModesAircraft
    scatter!(modeDampingRatiosAircraft[mode], modeFrequenciesAircraft[mode], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label=false)
    scatter!([modeDampingRatiosAircraft[mode][1]], [modeFrequenciesAircraft[mode][1]], c=colors[1], shape=mshape[1], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
for mode in 1:nModesWing
    scatter!(modeDampingRatiosWing[mode], modeFrequenciesWing[mode], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label=false)
    scatter!([modeDampingRatiosWing[mode][1]], [modeFrequenciesWing[mode][1]], c=colors[2], shape=mshape[2], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RLdampr)
savefig(string(absPath,"/cHALE_compare_stability_rootlocus_lambda_",λ,"_k2_",k2,".pdf"))

# Root locus - damping vs frequency
plt_RLdamp = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=(0.1,0.55))
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label=string("Aircraft - \$k_2=\$",k2))
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label=string("Trimmed wing- \$k_2=\$",k2))
for mode in 1:nModesAircraft
    scatter!(modeDampingsAircraft[mode], modeFrequenciesAircraft[mode], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label=false)
    scatter!([modeDampingsAircraft[mode][1]], [modeFrequenciesAircraft[mode][1]], c=colors[1], shape=mshape[1], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
for mode in 1:nModesWing
    scatter!(modeDampingsWing[mode], modeFrequenciesWing[mode], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label=false)
    scatter!([modeDampingsWing[mode][1]], [modeFrequenciesWing[mode][1]], c=colors[2], shape=mshape[2], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
if k2 == 0.015
    TIP1textPos = [0.75, 10]
    TIP2textPos = [0.75, 35]
    annotate!(TIP1textPos[1], TIP1textPos[2], text("1st T-IP", 12))
    annotate!(TIP2textPos[1], TIP2textPos[2], text("2nd T-IP", 12))
    quiver!([TIP1textPos[1]-0.5,TIP1textPos[1]-0.5,TIP2textPos[1]-0.25], [TIP1textPos[2]-2,TIP1textPos[2]+2,TIP2textPos[2]-2], quiver=([-0.65,-0.65,-0.55], [-2,4,-3]), arrow=:closed, linecolor=:black)
elseif k2 == 0.045
    TIP1textPos = [1, 8]
    TIP2textPos = [0.75, 35]
    annotate!(TIP1textPos[1], TIP1textPos[2], text("1st T-IP", 12))
    annotate!(TIP2textPos[1], TIP2textPos[2], text("2nd T-IP", 12))
    quiver!([TIP1textPos[1]-0.2,TIP1textPos[1],TIP2textPos[1]-0.25,TIP2textPos[1]-0.25], [TIP1textPos[2]+2,TIP1textPos[2]+2,TIP2textPos[2]-2,TIP2textPos[2]-2], quiver=([-0.65,0,-0.5,-0.3], [2,3,-2,-7]), arrow=:closed, linecolor=:black)
end
display(plt_RLdamp)
savefig(string(absPath,"/cHALE_compare_stability_rootlocus2_lambda_",λ,"_k2_",k2,".pdf"))

# V-g-f
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,120], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModesAircraft
    plot!(URange, modeFrequenciesAircraft[mode], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label=false)
end
for mode in 1:nModesWing
    plot!(URange, modeFrequenciesWing[mode], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.2], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModesAircraft
    plot!(URange, modeDampingsAircraft[mode]./modeFrequenciesAircraft[mode], c=colors[1], shape=mshape[1], ms=ms, msw=msw, label=false)
end
for mode in 1:nModesWing
    plot!(URange, modeDampingsWing[mode]./modeFrequenciesWing[mode], c=colors[2], shape=mshape[2], ms=ms, msw=msw, label=false)
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALE_compare_stability_Vgf_lambda_",λ,"_k2_",k2,".pdf"))

# Normalized deformed span at iminence of first flutter onset
plt_u3 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized OOP direction", xlims=[0,1], ylims=[-0.2,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(x1_0/L, x3_0/L, c=:gray, lw=lw, ls=:dot, label="Undeformed")
plot!(x1_def[flutterSpeedIndexAircraft]/L, x3_def[flutterSpeedIndexAircraft]/L, c=colors[1], lw=lw, label="At flutter - aircraft")
plot!(x1_def[flutterSpeedIndexWing]/L, x3_def[flutterSpeedIndexWing]/L, c=colors[2], lw=lw, label="At flutter - trimmed wing")
display(plt_u3)
savefig(string(absPath,"/cHALE_compare_stability_disp_lambda_",λ,"_k2_",k2,".pdf"))

# Flutter speeds
indA = 1:1:nModesAircraft
indW = 1:1:nModesWing
flutterSpeedAircraftAllSorted = sort(flutterSpeedAircraftAll)
flutterSpeedWingAllSorted = sort(flutterSpeedWingAll)
plt_Uf = plot(xlabel="Flutter occurrence", ylabel="Speed [m/s]", ylims=[0,URange[end]], tickfont=font(ts), guidefont=font(fs), xticks=indA)
plot!(indA[1:length(flutterSpeedAircraftAllSorted)],flutterSpeedAircraftAllSorted, c=colors[1], lw=lw, marker=(mshape[1], 5), msw=0, label=false)
plot!(indW[1:length(flutterSpeedWingAllSorted)],flutterSpeedWingAllSorted, c=colors[2], lw=lw, marker=(mshape[2], 5), label=false)
display(plt_Uf)
savefig(string(absPath,"/cHALE_compare_stability_flutterSpeeds_",λ,"_k2_",k2,".pdf"))

# Trim root angle of attack
plt_trimAoA = plot(xlabel="Airspeed [m/s]", ylabel="Wing root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=[-2,15], tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt_trimAoA)
savefig(string(absPath,"/cHALE_compare_stability_trimAoA_lambda_",λ,"_k2_",k2,".pdf"))

# Trim thrust
plt_trimThrust = plot(xlabel="Airspeed [m/s]", ylabel="Thrust [N]", xlims=[URange[1],URange[end]], ylims=[0,300], tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt_trimThrust)
savefig(string(absPath,"/cHALE_compare_stability_trimThrust_lambda_",λ,"_k2_",k2,".pdf"))

# Trim elevator deflection
plt_trimDelta = plot(xlabel="Airspeed [m/s]", ylabel="Elevator deflection [deg]", xlims=[URange[1],URange[end]], ylims=[-15,5], tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt_trimDelta)
savefig(string(absPath,"/cHALE_compare_stability_trimDelta_lambda_",λ,"_k2_",k2,".pdf"))

println("Finished cHALE_compare_stability.jl")