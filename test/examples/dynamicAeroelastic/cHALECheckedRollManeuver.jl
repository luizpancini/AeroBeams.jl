# # Coordinated turn maneuver of a HALE aircraft
# This example illustrates how to set up a dynamic analysis of an aircraft in free flight. For that we take a high-altitude long-endurance (HALE) aircraft model described by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738):

#md # ![](../assets/cHALE.png)
#md # *HALE model geometry* by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl).

# ### Problem setup
# Let's begin by setting up the variables of our problem.
using AeroBeams

## Altitude [m]
h = 20e3

## Airspeed [m/s]
U = 25

## Wing and stabilizers parasite drag and flag for induced drag
wingCd0 = stabsCd0 = 1e-2
hasInducedDrag = true

## Aerodynamic solver
aeroSolver = Indicial()

## Discretization
nElemWing = 20
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 4
#md nothing #hide

# ### Trim problem
# We have to first trim the aircraft at the specified flight condition. Our trim variables are the engine's thrust (modeled as a follower force at the root of the wing), and the elevator deflection. For the trim problem, we set a Newton-Raphson solver through [`create_NLsolver`](@ref create_NLsolver) with the pseudoinverse method for the linear step (`trimSolverMethod=:pinv`). We use the built-in function [`create_cHALE`](@ref create_cHALE) to streamline the model creation process. A dedicated [example](cHALEmodel.md) explains the inner workings of that function.
## System solver
NRtrim = create_NewtonRaphson(trimSolverMethod=:pinv)

## Model for trim problem
cHALEtrim,_ = create_cHALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,δElevIsTrimVariable=true,thrustIsTrimVariable=true)

## Create and solve trim problem
trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim)
solve!(trimProblem)

## Get wing, horizontal amd vertical tail elements
leftWingRootElem = div(nElemWing,2)
rightWingRootElem = leftWingRootElem + 1
leftHTRootElem = nElemWing + nElemTailBoom + div(nElemHorzStabilizer,2)
rightHTRootElem = leftHTRootElem + 1
leftVTRootElem = nElemWing + nElemTailBoom + nElemHorzStabilizer + div(nElemVertStabilizer,2)
rightVTRootElem = leftVTRootElem + 1

## Extract trim variables and outputs
trimSol = trimProblem.trimSolution
trimAeroVars = trimProblem.loadstepHistory.aeroVariables[end]
trimAoA = (trimAeroVars[leftWingRootElem].unsteadyKin.αₑ + trimAeroVars[rightWingRootElem].unsteadyKin.αₑ)/2
trimThrust = trimSol.trimLoads[1]
trimδ = trimSol.trimδ[1]
println("Trim variables/outputs: AoA = $(round(rad2deg(trimAoA), digits=2)) deg, T = $(round(trimThrust, digits=2)), δ = $(round(rad2deg(trimδ), digits=2)) deg") #src
#md nothing #hide

# ### Dynamic problem
# The coordinate turn maneuver investigated is defined by checked aileron and rudder deflections, linearly ramped. The time-dependent aileron and rudder profiles, `δAil` and `δRudd`, are passed as arguments to create the model for the dynamic problem.
#md using Suppressor #hide
## Set checked aileron and rudder deflection profiles
Δδ = deg2rad(-15)
tδinit = 1
tδramp = 5
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp
δRudd = δAil = t -> ifelse(
    t <= tδinit, 
    0,
    ifelse(
        t <= tδpeak, 
        0 + Δδ * ((t-tδinit) / (tδpeak-tδinit)),
        ifelse(
            t <= tδfinal, 
            0 + Δδ - Δδ * ((t-tδpeak) / (tδfinal-tδpeak)),
            0
        )
    )
)

## Model for dynamic problem
cHALEdyn,_ = create_cHALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,δElev=trimδ,δAil=δAil,δRudd=δRudd,thrust=trimThrust)
#md nothing #hide

## Time variables
Δt = 5e-2
tf = 30

## Set system solver for dynamic problem
NRdyn = create_NewtonRaphson()

## Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=cHALEdyn,x0=trimProblem.systemArrays.x[cHALEtrim.assemblyData.rangeNominal],finalTime=tf,Δt=Δt,systemSolver=NRdyn)
#md @suppress begin #hide
solve!(dynamicProblem)
#md end #hide
#md nothing #hide

# ### Post-processing
# The first post-processing step is to unpack the desired outputs.
## Unpack numerical solution
sol = dynamicProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
wingAoA = rad2deg.([(aeroVars[i][leftWingRootElem].unsteadyKin.αₑ + aeroVars[i][rightWingRootElem].unsteadyKin.αₑ)/2 for i in eachindex(t)])
htAoA = rad2deg.([(aeroVars[i][leftHTRootElem].unsteadyKin.αₑ + aeroVars[i][rightHTRootElem].unsteadyKin.αₑ)/2 for i in eachindex(t)])
vtAoA = rad2deg.([(aeroVars[i][leftVTRootElem].unsteadyKin.αₑ + aeroVars[i][rightVTRootElem].unsteadyKin.αₑ)/2 for i in eachindex(t)])
leftWingtipAoA = rad2deg.([aeroVars[i][1].unsteadyKin.αₑ  for i in eachindex(t)])
rightWingtipAoA = rad2deg.([aeroVars[i][nElemWing].unsteadyKin.αₑ for i in eachindex(t)])
Δu1 = [(cHALEdyn.kinematicsBasisA.u_AVec[i][1] + states[i][leftWingRootElem].nodal.u2[1]) for i in eachindex(t)] .- states[1][leftWingRootElem].nodal.u2[1]
Δu2 = [(cHALEdyn.kinematicsBasisA.u_AVec[i][2] + states[i][leftWingRootElem].nodal.u2[2]) for i in eachindex(t)] .- states[1][leftWingRootElem].nodal.u2[2]
Δu3 = [(cHALEdyn.kinematicsBasisA.u_AVec[i][3] + states[i][leftWingRootElem].nodal.u2[3]) for i in eachindex(t)] .- states[1][leftWingRootElem].nodal.u2[3]
airspeed = [norm(aeroVars[i][leftWingRootElem].steadyKin.U + aeroVars[i][rightWingRootElem].steadyKin.U)/2 for i in eachindex(t)]
#md nothing #hide

#md # We are now ready to visualize the results through the changes in lateral distance (``x_1``), forward distance (``x_2``) and altitude (``x_3``), angles of attack and airspeed.
#md gr()
#md ENV["GKSwstype"] = "100" #hide

## Lateral (x1) distance
#md plt1 = plot(xlabel="Time [s]", ylabel="Lateral distance [m]")
#md plot!(t, Δu1, c=:black, lw=2, label=false)
#md savefig("cHALECheckedRollManeuver_ldist.svg") #hide

## Forward (x2) distance
#md plt2 = plot(xlabel="Time [s]", ylabel="Forward distance [m]")
#md plot!(t, Δu2, c=:black, lw=2, label=false)
#md savefig("cHALECheckedRollManeuver_fdist.svg") #hide

## Altitude (Δx3)
#md plt3 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
#md plot!(t, Δu3, c=:black, lw=2, label=false)
#md savefig("cHALECheckedRollManeuver_altitude.svg") #hide

## AoA
#md plt4 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
#md plot!(t, wingAoA*180/π, c=:black, lw=2, label="Wing (root)")
#md plot!(t, leftWingtipAoA*180/π, c=:green, lw=2, label="Wing (left tip)")
#md plot!(t, rightWingtipAoA*180/π, c=:red, lw=2, label="Wing (right tip)")
#md plot!(t, htAoA*180/π, c=:blue, lw=2, label="HT (root)")
#md plot!(t, vtAoA*180/π, c=:cyan, lw=2, label="VT (root)")
#md savefig("cHALECheckedRollManeuver_AoA.svg") #hide

## Airspeed
#md plt5 = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]")
#md plot!(t, airspeed, c=:black, lw=2, label=false)
#md savefig("cHALECheckedRollManeuver_airspeed.svg") #hide
#md nothing #hide

#md # ![](cHALECheckedRollManeuver_ldist.svg)
#md # ![](cHALECheckedRollManeuver_fdist.svg)
#md # ![](cHALECheckedRollManeuver_altitude.svg)
#md # ![](cHALECheckedRollManeuver_AoA.svg)
#md # ![](cHALECheckedRollManeuver_airspeed.svg)

#md # Let's also visualize the motion of the aircraft through the view of an observer that follows the aircraft, keeping it at the center of the plot box. For that, we use the [`plot_dynamic_deformation`](@ref plot_dynamic_deformation) function with `follow=:assembly` and other appropriate arguments.
#md @suppress begin #hide
#md plot_dynamic_deformation(dynamicProblem,follow=:assembly,view=(30,15),plotDistLoads=false,plotFrequency=4,fps=24,plotLimits=([-20,20],[-20,20],[-20,20]),save=true,savePath="/docs/build/literate/cHALECheckedRollManeuver_motion.gif")
#md nothing #hide
#md end #hide
#md # ![](cHALECheckedRollManeuver_motion.gif)
#md nothing #hide

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALECheckedRollManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:assembly,view=(30,15),plotDistLoads=false,plotFrequency=4,fps=24,plotLimits=([-20,20],[-20,20],[-20,20]),save=true,savePath=string(relPath,"/cHALECheckedRollManeuver_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Lateral distance
plt1 = plot(xlabel="Time [s]", ylabel="Lateral distance [m]")
plot!(t, Δu1, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALECheckedRollManeuver_ldist.pdf"))

# Forward distance
plt1 = plot(xlabel="Time [s]", ylabel="Forward distance [m]")
plot!(t, Δu2, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALECheckedRollManeuver_fdist.pdf"))

# Altitude
plt2 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALECheckedRollManeuver_altitude.pdf"))

# AoA
plt3 = plot(xlabel="Time [s]", ylabel="Angle of attack [deg]")
plot!(t, wingAoA, c=:black, lw=lw, label="Wing (center)")
plot!(t, leftWingtipAoA, c=:green, lw=lw, label="Wing (left tip)")
plot!(t, rightWingtipAoA, c=:red, lw=lw, label="Wing (right tip.)")
plot!(t, htAoA, c=:blue, lw=lw, label="HT (center)")
plot!(t, vtAoA, c=:cyan, lw=lw, label="VT (root)")
display(plt3)
savefig(string(absPath,"/cHALECheckedRollManeuver_AoA.pdf"))

# Airspeed
plt4 = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]")
plot!(t, airspeed, c=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/cHALECheckedRollManeuver_airspeed.pdf"))

println("Finished cHALECheckedRollManeuver.jl") #src