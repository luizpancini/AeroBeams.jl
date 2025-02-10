using AeroBeams

# Aerodynamic solvers
aeroSolvers = [BLi()]

# Wing airfoil
wingAirfoil = NACA23012A

# Option for reduced chord
reducedChord = true

# TF to include beam pods and number of elements
beamPods = false

# Option to set payload on wing
payloadOnWing = true

# Number of elements of beams
nElemStraightSemispan = 10
nElemDihedralSemispan = 5
nElemPod = 2

# Payload [lb]
P = 0

# Airspeed [m/s]
U = 40*0.3048

# Elevator profile variables
Δδ = 8*π/180
tδinit = 1
tδramp = 1
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp

# Time variables
Δt = [5e-2, 1e-2]
tf = 20

# Set NR system solver for trim problem
relaxFactor = 0.5
maxIter = 50
rtol = 1e-12
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,relativeTolerance=rtol,displayStatus=false)

# Set NR system solver for dynamic problem
maxIter = 100
alwaysUpdateJacobian = false
minConvRateAeroJacUpdate = 1.2
minConvRateJacUpdate = 1.2
rtol = 1e-8
NRdyn = create_NewtonRaphson(maximumIterations=maxIter,alwaysUpdateJacobian=alwaysUpdateJacobian,minConvRateAeroJacUpdate=minConvRateAeroJacUpdate,minConvRateJacUpdate=minConvRateJacUpdate,relativeTolerance=rtol,displayStatus=false)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(aeroSolvers))
dynamicProblem = Array{DynamicProblem}(undef,length(aeroSolvers))
t = Array{Vector{Float64}}(undef,length(aeroSolvers))
rootAoA = Array{Vector{Float64}}(undef,length(aeroSolvers))
Δu3 = Array{Vector{Float64}}(undef,length(aeroSolvers))
airspeed = Array{Vector{Float64}}(undef,length(aeroSolvers))
root_κ2 = Array{Vector{Float64}}(undef,length(aeroSolvers))
χ = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
χdot = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Model for trim problem
    heliosTrim,_ = create_Helios(reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,aeroSolver=aeroSolver,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
    # Create and solve trim problem
    trimProblem[i] = create_TrimProblem(model=heliosTrim,systemSolver=NRtrim)
    solve!(trimProblem[i])
    # Extract trim variables
    trimThrust = trimProblem[i].x[end-1]*trimProblem[i].model.forceScaling
    trimδ = trimProblem[i].x[end]
    println("Trim variables: T = $(trimThrust), δ = $(trimδ*180/π)")
    # Set checked elevator deflection profile
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
    heliosDynamic,midSpanElem,_ = create_Helios(reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,aeroSolver=aeroSolver,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,payloadPounds=P,airspeed=U,δ=δ,thrust=trimThrust)
    # Create and solve dynamic problem
    dynamicProblem[i] = create_DynamicProblem(model=heliosDynamic,x0=trimProblem[i].x[1:end-2],finalTime=tf,Δt=Δt[i],skipInitialStatesUpdate=true,systemSolver=NRdyn)
    solve!(dynamicProblem[i])
    # Outputs
    t[i] = dynamicProblem[i].savedTimeVector
    Nt = length(t[i])
    rootAoA[i] = [dynamicProblem[i].aeroVariablesOverTime[j][midSpanElem].flowAnglesAndRates.αₑ for j in 1:Nt]
    Δu3[i] = [dynamicProblem[i].nodalStatesOverTime[j][midSpanElem].u_n2[3] for j in 1:Nt] .- dynamicProblem[i].nodalStatesOverTime[1][midSpanElem].u_n2[3]
    airspeed[i] = [dynamicProblem[i].aeroVariablesOverTime[j][midSpanElem].flowVelocitiesAndRates.U∞ for j in 1:Nt]
    root_κ2[i] = [dynamicProblem[i].compElementalStatesOverTime[j][midSpanElem].κ[2] for j in 1:Nt]
    χ[i] = [dynamicProblem[i].elementalStatesOverTime[j][midSpanElem].χ for j in 1:Nt]
    χdot[i] = [dynamicProblem[i].elementalStatesRatesOverTime[j][midSpanElem].χdot for j in 1:Nt]
end

# Set paths
relPath = "/dev/helios/figures/heliosFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
solversLabels = ["Linear", "Dynamic stall"]
solversLS = [:dash, :solid]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)
stateColors = cgrad(:rainbow, 8, categorical=true)

# Root angle of attack
plt_AoA = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]", tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=12)
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(t[i], rootAoA[i]*180/π, c=colors[i], lw=lw, ls=solversLS[i], label=solversLabels[i])
end
display(plt_AoA)
savefig(string(absPath,"/heliosPitchManeuver_AoA.pdf"))

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(t[i], airspeed[i], c=colors[i], lw=lw, ls=solversLS[i], label=false)
end
display(plt_U)
savefig(string(absPath,"/heliosPitchManeuver_U.pdf"))

# Elevation change
plt_deltaz = plot(xlabel="Time [s]", ylabel="Elevation change [m]", tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(t[i], Δu3[i], c=colors[i], lw=lw, ls=solversLS[i], label=false)
end
display(plt_deltaz)
savefig(string(absPath,"/heliosPitchManeuver_deltaz.pdf"))

# Root bending curvature
plt_k2 = plot(xlabel="Time [s]", ylabel="Root bending curvature [1/m]", tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(t[i], root_κ2[i], c=colors[i], lw=lw, ls=solversLS[i], label=false)
end
display(plt_k2)
savefig(string(absPath,"/heliosPitchManeuver_k2.pdf"))

# States
plt_χ = plot(xlabel="Time [s]", ylabel="States", tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for j in eachindex(χ[i][1])
        χ_now = [χ[i][k][j] for k in eachindex(t[i])]
        plot!(t[i], χ_now, c=stateColors[j], lw=lw, label=false)
    end
end
display(plt_χ)

# States' rates
plt_χdot = plot(xlabel="Time [s]", ylabel="States rates", tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for j in eachindex(χdot[i][1])
        χdot_now = [χdot[i][k][j] for k in eachindex(t[i])]
        plot!(t[i], χdot_now, c=stateColors[j], lw=lw, label=false)
    end
end
display(plt_χdot)

println("Finished heliosPitchManeuver.jl")