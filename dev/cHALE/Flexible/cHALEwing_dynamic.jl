using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Altitude [m]
h = 20e3

# Airspeed [m/s]
U = 36

# Stiffness factor
λ = 1

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Bending pre-curvature
k2 = 0.045

# Discretization
if λ == 1
    nElemWing = 80
elseif λ > 1
    nElemWing = 40
end
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# Number of modes for eigenproblem
nModes = 20

# Flag to solve preliminary trim problem at smaller airspeed
solvePrelimTrim = true
UprelimTrim = [20]

# System solver for trim and eigen problems
relaxFactor = 0.5
maxIter = 100
σ0t = 1
σ0e = 1
relTol = 1e-8
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,desiredIterations=div(maxIter,4),initialLoadFactor=σ0t,relativeTolerance=relTol,displayStatus=true)
NReigen = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,desiredIterations=div(maxIter,4),initialLoadFactor=σ0e,relativeTolerance=relTol,displayStatus=true)

# Model for trim problem
cHALEtrim,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=true,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)

# Solve trim problem at smaller airspeed for better initial guess, if applicable
global x0Trim = zeros(0)
if solvePrelimTrim
    cHALEtrim.skipValidationMotionBasisA = false
    for Uprelim in UprelimTrim
        println("Solving preliminary trim problem at U=$Uprelim")
        set_motion_basis_A!(model=cHALEtrim,v_A=[0;Uprelim;0])
        prelimTrimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(prelimTrimProblem)
        global x0Trim = prelimTrimProblem.x
    end
    set_motion_basis_A!(model=cHALEtrim,v_A=[0;U;0])
end

# Create and solve trim problem
trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
println("Solving trim problem")
solve!(trimProblem)

# Extract trim variables and outputs
trimAoA = (trimProblem.aeroVariablesOverσ[end][cHALEtrim.beams[1].elementRange[end]].flowAnglesAndRates.αₑ + trimProblem.aeroVariablesOverσ[end][cHALEtrim.beams[2].elementRange[1]].flowAnglesAndRates.αₑ)/2
trimThrust = trimProblem.x[end-1]*trimProblem.model.forceScaling
trimδ = trimProblem.x[end]
println("Trim outputs: AoA = $(trimAoA*180/π), T = $(trimThrust), δ = $(trimδ*180/π)")

# Tip disturbance force
δF = 50
t₀ = 1
tipF3 = t -> ifelse(
        t <= t₀,
        δF*t/t₀,
        ifelse(
            t <= 2*t₀,
            δF*t₀ - δF*(t-t₀),
            0
        )
)

# Model for eigen problem
wingEigen,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA)

# Create and solve eigen problem
eigenProblem = create_EigenProblem(model=wingEigen,nModes=nModes,frequencyFilterLimits=[0,Inf],systemSolver=NReigen)
println("Solving eigenproblem")
solve!(eigenProblem)

# Frequencies and dampings
freqs = eigenProblem.frequenciesOscillatory
damps = round_off!(eigenProblem.dampingsOscillatory,1e-8)

# Show whether flutter is expected from eigenanalysis
if any(x->x>0,damps)
    ind = findall(x->x>0,damps)
    for i in ind
        println("Flutter is expected, damps[",i,"] = $(damps[i]), ωf = $(freqs[i])")
    end
else
    println("Flutter is NOT expected")
end

# Model for dynamic problem
wingDynamic,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA,tipF3=tipF3)

# Time variables
Δt = 5e-2
tf = 120

# Set NR system solver for dynamic problem
maxIter = 50
NRdyn = create_NewtonRaphson(maximumIterations=maxIter)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=wingDynamic,finalTime=tf,Δt=Δt,systemSolver=NRdyn,skipInitialStatesUpdate=true,x0=eigenProblem.x)
println("Solving dynamic problem")
solve!(dynamicProblem)

# Unpack numerical solution
t = dynamicProblem.savedTimeVector
tipAoA = [dynamicProblem.aeroVariablesOverTime[i][end].flowAnglesAndRates.αₑ for i in 1:length(t)]

# Set paths
relPathFig = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_dynamic"
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALEwing_dynamic"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
mkpath(absPathFig)
mkpath(absPathData)

# Plots
using Plots, DelimitedFiles
lw = 2
gr()

# Tip AoA
plt_AoA = plot(xlabel="Time [s]", ylabel="Normalized tip angle of attack", xticks=0:30:tf)
plot!(t, tipAoA./tipAoA[1], c=:black, lw=lw, label=false)
display(plt_AoA)
savefig(string(absPathFig,string("/cHALEwing_dynamic_AoA_lambda",λ,"_U",U,"_k2",k2,".pdf")))

# Save arrays
writedlm(string(absPathData,"/cHALEwing_dynamic_t_lambda",λ,"_U",U,"_k2",k2,".txt"), t)
writedlm(string(absPathData,"/cHALEwing_dynamic_AoA_lambda",λ,"_U",U,"_k2",k2,".txt"), tipAoA)

# Estimated oscillation frequencies from consecutive local maxima
local_maxima_indices = findall(x -> x < 0, diff(sign.(diff(tipAoA))))
t_local_maxima_indices = t[local_maxima_indices]
freqs_from_maxima = []
for i in 2:length(t_local_maxima_indices)
    push!(freqs_from_maxima,2π/(t_local_maxima_indices[i]-t_local_maxima_indices[i-1]))
end

# # Show interactive plot
# plotlyjs()
# plt_AoA2 = plot(xlabel="Time [s]", ylabel="Normalized tip angle of attack", xticks=0:30:tf)
# plot!(t, tipAoA./tipAoA[1], c=:black, lw=lw, label=false)
# display(plt_AoA2)

println("Finished cHALEwing_dynamic.jl")