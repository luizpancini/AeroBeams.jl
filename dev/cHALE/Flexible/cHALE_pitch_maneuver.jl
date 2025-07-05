using AeroBeams, JLD2

# Aerodynamic solver
aeroSolver = Indicial()

# Altitude [m]
h = 20e3

# Airspeed [m/s]
U = 45

# Stiffness factor
λ = 2

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Flag to solve preliminary trim problem at smaller airspeed
solvePrelimTrim = true
UprelimTrim = [20,30]

# Bending pre-curvature
k2 = 0.030

# Discretization
if λ == 1
    nElemWing = 80
elseif λ > 1
    nElemWing = 40
end
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# Number of modes for eigenproblem
nModes = 25

# Attachment springs for eigenproblem
μu = μp = 1e-2
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=μu*[1; 1; 1],kp=μp*[1; 1; 1])

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1
relTol = 1e-8
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,relativeTolerance=relTol,displayStatus=true)

# Model for trim problem without springs
cHALEtrim,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)

# Solve trim problem at smaller airspeed for better initial guess, if applicable
global x0 = zeros(0)
if solvePrelimTrim
    cHALEtrim.skipValidationMotionBasisA = false
    for Uprelim in UprelimTrim
        println("Solving preliminary trim problem at U=$Uprelim")
        set_motion_basis_A!(model=cHALEtrim,v_A=[0;Uprelim;0])
        prelimTrimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR,x0=x0)
        solve!(prelimTrimProblem)
        global x0 = prelimTrimProblem.x
    end
    set_motion_basis_A!(model=cHALEtrim,v_A=[0;U;0])
end

# Create and solve trim problem
trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR,x0=x0)
println("Solving trim problem")
solve!(trimProblem)

# Extract trim variables and outputs
trimAoA = (trimProblem.aeroVariablesOverσ[end][cHALEtrim.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblem.aeroVariablesOverσ[end][cHALEtrim.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]
println("Trim outputs: AoA = $(trimAoA*180/π), T = $(trimThrust), δ = $(trimδ*180/π)")

# Model for trim problem with springs
cHALEtrimSpringed,_,_,tailBoomSpringed,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)

# Add springs
add_springs_to_beam!(beam=tailBoomSpringed,springs=[spring1])

# Update model
cHALEtrimSpringed.skipValidationMotionBasisA = true
update_model!(cHALEtrimSpringed)

# Create and solve trim problem with springs
trimProblemSpringed = create_TrimProblem(model=cHALEtrimSpringed,systemSolver=NR,x0=trimProblem.x)
println("Solving trim problem with springs")
solve!(trimProblemSpringed)

# Compare trim outputs
trimAoASpringed = (trimProblemSpringed.aeroVariablesOverσ[end][cHALEtrimSpringed.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblemSpringed.aeroVariablesOverσ[end][cHALEtrimSpringed.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
trimThrustSpringed = trimProblemSpringed.x[end-1]*trimProblemSpringed.model.forceScaling
trimδSpringed = trimProblemSpringed.x[end]
println("Trim outputs ratios springed/nominal: AoA = $(trimAoASpringed/trimAoA), T = $(trimThrustSpringed/trimThrust), δ = $(trimδSpringed/trimδ)")

# Model for eigen problem (with springs)
cHALEeigen,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ,thrust=trimThrust,k2=k2,hasInducedDrag=hasInducedDrag)

# Create and solve eigen problem
println("Solving eigenproblem")
eigenProblem = create_EigenProblem(model=cHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64],refTrimProblem=trimProblemSpringed)
solve_eigen!(eigenProblem)

# Frequencies and dampings
freqs = eigenProblem.frequenciesOscillatory
damps = round_off!(eigenProblem.dampingsOscillatory,1e-8)
dampRatio = damps./freqs

# Show whether flutter is expected from eigenanalysis
if any(x->x>0,damps)
    ind = findall(x->x>0,damps)
    for i in ind
        println("Flutter is expected, damps[",i,"] = $(damps[i]), ωf = $(freqs[i])")
    end
else
    println("Flutter is NOT expected")
end

# Set checked elevator deflection profile for dynamic problem
Δδ = -.1*π/180
tδinit = 0.5
tδramp = 0.5
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp
δ = t -> ifelse(
    t <= tδinit, 
    trimδ,
    ifelse(
        t <= tδpeak, 
        trimδ + Δδ * ((t-tδinit) / (tδpeak-tδinit)),
        ifelse(
            t <= tδfinal, 
            trimδ + Δδ - Δδ * ((t-tδpeak) / (tδfinal-tδpeak)),
            trimδ
        )
    )
)

# Model for dynamic problem
cHALEdynamic,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,k2=k2,δElev=δ,thrust=trimThrust)

# Time variables
Δt = 5e-2
tf = 120

# Set NR system solver for dynamic problem
maxIter = 50
NRdyn = create_NewtonRaphson(maximumIterations=maxIter)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=cHALEdynamic,x0=trimProblem.x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true,systemSolver=NRdyn)
solve!(dynamicProblem)

# Get wing root elements
lRootElem = div(nElemWing,2)
rRootElem = lRootElem+1

# Unpack numerical solution
t = dynamicProblem.savedTimeVector
wingAoA = [(dynamicProblem.aeroVariablesOverTime[i][lRootElem].flowAnglesAndRates.αₑ+dynamicProblem.aeroVariablesOverTime[i][rRootElem].flowAnglesAndRates.αₑ)/2 for i in 1:length(t)]
airspeed = [(dynamicProblem.aeroVariablesOverTime[i][lRootElem].flowVelocitiesAndRates.U∞+dynamicProblem.aeroVariablesOverTime[i][rRootElem].flowVelocitiesAndRates.U∞)/2 for i in 1:length(t)]

# Set paths
relPathFig = "/dev/cHALE/Flexible/outputs/figures/cHALE_pitch_maneuver"
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALE_pitch_maneuver"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
mkpath(absPathFig)
mkpath(absPathData)

# Plots
using Plots, DelimitedFiles
lw = 2
gr()

# Wing root AoA
plt_AoA = plot(xlabel="Time [s]", ylabel="Normalized wing root angle of attack", xticks=0:30:tf)
plot!(t, wingAoA./wingAoA[1], c=:black, lw=lw, label=false)
display(plt_AoA)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_AoA_lambda",λ,"_U",U,"_k2",k2,".pdf")))

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Normalized airspeed", xticks=0:30:tf)
plot!(t, airspeed./airspeed[1], c=:black, lw=lw, label=false)
display(plt_U)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_airspeed_lambda",λ,"_U",U,"_k2",k2,".pdf")))

# Save arrays
writedlm(string(absPathData,"/cHALE_pitch_maneuver_t_lambda",λ,"_U",U,"_k2",k2,".txt"), t)
writedlm(string(absPathData,"/cHALE_pitch_maneuver_AoA_lambda",λ,"_U",U,"_k2",k2,".txt"), wingAoA)

# Estimated oscillation frequencies from consecutive local maxima
local_maxima_indices = findall(x -> x < 0, diff(sign.(diff(wingAoA))))
t_local_maxima_indices = t[local_maxima_indices]
freqs_from_maxima = []
for i in 2:length(t_local_maxima_indices)
    push!(freqs_from_maxima,2π/(t_local_maxima_indices[i]-t_local_maxima_indices[i-1]))
end

# # Show interactive plot
# plotlyjs()
# plt_AoA2 = plot(xlabel="Time [s]", ylabel="Normalized wing root angle of attack", xticks=0:30:tf)
# plot!(t, wingAoA./wingAoA[1], c=:black, lw=lw, label=false)
# display(plt_AoA2)

println("Finished cHALE_pitch_maneuver.jl")