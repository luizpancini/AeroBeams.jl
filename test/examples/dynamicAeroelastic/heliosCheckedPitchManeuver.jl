using AeroBeams

# Stiffness factor range
λRange = [1,2,10]

# Payload [lb]
P = 200

# Airspeed [m/s]
U = feet2meters(40)

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = false

# TF to include beam pods and number of elements
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Number of elements of beams
nElemStraightSemispan = 10
nElemDihedralSemispan = 5
nElemPod = 2

# Aerodynamic solver
aeroSolver = Indicial()

# Time variables
Δt = 5e-2
tf = 60

# Checked elevator deflection profile
Δδ = deg2rad(-5)
tδinit = 5
tδramp = 2
tδpeak = tδinit + tδramp
tδfinal = tδpeak + tδramp

# System solver for trim problem
NRtrim = create_NewtonRaphson(trimSolverMethod=:pinv,absTol=1e-7,displayStatus=true)

# System solver for dynamic problem
NRdyn = create_NLsolver()

# Initialize outputs
dynamicProblem = Array{DynamicProblem}(undef, length(λRange))
trimThrust = Array{Float64}(undef, length(λRange))
trimδ = Array{Float64}(undef, length(λRange))
t = Array{Vector{Float64}}(undef, length(λRange))
rootAoA = Array{Vector{Float64}}(undef, length(λRange))
Δu3 = Array{Vector{Float64}}(undef, length(λRange))
airspeed = Array{Vector{Float64}}(undef, length(λRange))

# Loop stiffness factor range
for (i, λ) in enumerate(λRange)
    println("Running for λ = $λ")
    # Model for trim problem
    heliosTrim,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing)
    # Create and solve trim problem
    trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NRtrim)
    solve!(trimProblem)
    # Extract trim variables
    trimSol = trimProblem.trimSolution
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = trimSol.trimδ[1]
    println("Trim variables: T = $(round(trimThrust[i], digits=2)), δ = $(round(rad2deg(trimδ[i]), digits=2)) deg")
    # Set checked elevator deflection profile
    δ = t -> ifelse(
        t <= tδinit, 
        trimδ[i],
        ifelse(
            t <= tδpeak, 
            trimδ[i] + Δδ * ((t-tδinit) / (tδpeak-tδinit)),
            ifelse(
                t <= tδfinal, 
                trimδ[i] + Δδ - Δδ * ((t-tδpeak) / (tδfinal-tδpeak)),
                trimδ[i]
            )
        )
    )
    # Model for dynamic problem
    heliosDynamic,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=δ,thrust=trimThrust[i],reducedChord=reducedChord,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing)
    # Create and solve dynamic problem
    dynamicProblem[i] = create_DynamicProblem(model=heliosDynamic,x0=trimProblem.systemArrays.x[heliosTrim.assemblyData.rangeNominal],finalTime=tf,Δt=Δt,systemSolver=NRdyn)
    solve!(dynamicProblem[i])
    # Unpack numerical solution
    sol = dynamicProblem[i].solutionHistory
    t[i] = sol.t
    states = sol.states
    aeroVars = sol.aeroVariables
    rootAoA[i] = rad2deg.([aeroVars[j][midSpanElem].unsteadyKin.αₑ for j in eachindex(t[i])])
    Δu3[i] = [states[j][midSpanElem].nodal.u2[3] for j in eachindex(t[i])] .- states[1][midSpanElem].nodal.u2[3]
    airspeed[i] = [norm(aeroVars[j][midSpanElem].steadyKin.U) for j in eachindex(t[i])]
end

## Plots
# Set paths
relPath = "/test/outputs/figures/heliosCheckedPitchManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(λRange), categorical=true)
lw = 2
gr()

# Animation
anim = plot_dynamic_deformations(dynamicProblem,legendEntries=["\$\\lambda=$(λ)\$" for λ in λRange],view=(30,15),follow=:assembly,plotDistLoads=false,plotFrequency=10,plotLimits=([-45,45],[-45,45],[-45,45]),save=true,savePath=string(relPath,"/heliosCheckedPitchManeuver_deformation.gif"),displayProgress=true)
display(anim)

# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
for (i, λ) in enumerate(λRange)
    plot!(t[i], Δu3[i], color=colors[i], lw=lw, label="\$\\lambda=$(λ)\$")
end
display(plt1)
savefig(string(absPath,"/heliosCheckedPitchManeuver_altitude.pdf"))

# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
for (i, λ) in enumerate(λRange)
    plot!(t[i], rootAoA[i], color=colors[i], lw=lw, label="\$\\lambda=$(λ)\$")
end
display(plt2)
savefig(string(absPath,"/heliosCheckedPitchManeuver_rootAoA.pdf"))

println("Finished heliosCheckedPitchManeuver.jl")