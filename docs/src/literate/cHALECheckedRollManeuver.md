```@meta
EditURL = "../../../test/examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl"
```

# Coordinated turn maneuver of a HALE aircraft
This example illustrates how to set up a dynamic analysis of an aircraft in free flight. For that we take a high-altitude long-endurance (HALE) aircraft model described by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738):

![](../assets/cHALE.png)
*HALE model geometry* by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738)

!!! tip
    The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl).

### Problem setup
Let's begin by setting up the variables of our problem.

````@example cHALECheckedRollManeuver
using AeroBeams

# Altitude [m]
h = 20e3

# Airspeed [m/s]
U = 25

# Wing and stabilizers parasite drag and flag for induced drag
wingCd0 = stabsCd0 = 1e-2
hasInducedDrag = true

# Aerodynamic solver
aeroSolver = Indicial()

# Discretization
nElemWing = 20
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 4
nothing #hide
````

### Trim problem
We have to first trim the aircraft at the specified flight condition. Our trim variables are the engine's thrust (modeled as a follower force at the root of the wing), and the elevator deflection. For the trim problem, we set a Newton-Raphson solver through [`create_NLsolver`](@ref create_NLsolver) with the pseudoinverse method for the linear step (`trimSolverMethod=:pinv`). We use the built-in function [`create_cHALE`](@ref create_cHALE) to streamline the model creation process. A dedicated [example](cHALEmodel.md) explains the inner workings of that function.

````@example cHALECheckedRollManeuver
# System solver
NRtrim = create_NewtonRaphson(trimSolverMethod=:pinv)

# Model for trim problem
cHALEtrim,_ = create_cHALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,δElevIsTrimVariable=true,thrustIsTrimVariable=true)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim)
solve!(trimProblem)

# Get wing, horizontal amd vertical tail elements
leftWingRootElem = div(nElemWing,2)
rightWingRootElem = leftWingRootElem + 1
leftHTRootElem = nElemWing + nElemTailBoom + div(nElemHorzStabilizer,2)
rightHTRootElem = leftHTRootElem + 1
leftVTRootElem = nElemWing + nElemTailBoom + nElemHorzStabilizer + div(nElemVertStabilizer,2)
rightVTRootElem = leftVTRootElem + 1

# Extract trim variables and outputs
trimSol = trimProblem.trimSolution
trimAeroVars = trimProblem.loadstepHistory.aeroVariables[end]
trimAoA = (trimAeroVars[leftWingRootElem].unsteadyKin.αₑ + trimAeroVars[rightWingRootElem].unsteadyKin.αₑ)/2
trimThrust = trimSol.trimLoads[1]
trimδ = trimSol.trimδ[1]
nothing #hide
````

### Dynamic problem
The coordinate turn maneuver investigated is defined by checked aileron and rudder deflections, linearly ramped. The time-dependent aileron and rudder profiles, `δAil` and `δRudd`, are passed as arguments to create the model for the dynamic problem.

````@example cHALECheckedRollManeuver
using Suppressor #hide
# Set checked aileron and rudder deflection profiles
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

# Model for dynamic problem
cHALEdyn,_ = create_cHALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,δElev=trimδ,δAil=δAil,δRudd=δRudd,thrust=trimThrust)
nothing #hide

# Time variables
Δt = 5e-2
tf = 30

# Set system solver for dynamic problem
NRdyn = create_NewtonRaphson()

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=cHALEdyn,x0=trimProblem.systemArrays.x[cHALEtrim.assemblyData.rangeNominal],finalTime=tf,Δt=Δt,systemSolver=NRdyn)
@suppress begin #hide
solve!(dynamicProblem)
end #hide
nothing #hide
````

### Post-processing
The first post-processing step is to unpack the desired outputs.

````@example cHALECheckedRollManeuver
# Unpack numerical solution
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
nothing #hide
````

We are now ready to visualize the results through the changes in lateral distance (``x_1``), forward distance (``x_2``) and altitude (``x_3``), angles of attack and airspeed.

````@example cHALECheckedRollManeuver
gr()
ENV["GKSwstype"] = "100" #hide

# Lateral (x1) distance
plt1 = plot(xlabel="Time [s]", ylabel="Lateral distance [m]")
plot!(t, Δu1, c=:black, lw=2, label=false)
savefig("cHALECheckedRollManeuver_ldist.svg") #hide

# Forward (x2) distance
plt2 = plot(xlabel="Time [s]", ylabel="Forward distance [m]")
plot!(t, Δu2, c=:black, lw=2, label=false)
savefig("cHALECheckedRollManeuver_fdist.svg") #hide

# Altitude (Δx3)
plt3 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, c=:black, lw=2, label=false)
savefig("cHALECheckedRollManeuver_altitude.svg") #hide

# AoA
plt4 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, wingAoA*180/π, c=:black, lw=2, label="Wing (root)")
plot!(t, leftWingtipAoA*180/π, c=:green, lw=2, label="Wing (left tip)")
plot!(t, rightWingtipAoA*180/π, c=:red, lw=2, label="Wing (right tip)")
plot!(t, htAoA*180/π, c=:blue, lw=2, label="HT (root)")
plot!(t, vtAoA*180/π, c=:cyan, lw=2, label="VT (root)")
savefig("cHALECheckedRollManeuver_AoA.svg") #hide

# Airspeed
plt5 = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]")
plot!(t, airspeed, c=:black, lw=2, label=false)
savefig("cHALECheckedRollManeuver_airspeed.svg") #hide
nothing #hide
````

![](cHALECheckedRollManeuver_ldist.svg)
![](cHALECheckedRollManeuver_fdist.svg)
![](cHALECheckedRollManeuver_altitude.svg)
![](cHALECheckedRollManeuver_AoA.svg)
![](cHALECheckedRollManeuver_airspeed.svg)

Let's also visualize the motion of the aircraft through the view of an observer that follows the aircraft, keeping it at the center of the plot box. For that, we use the [`plot_dynamic_deformation`](@ref plot_dynamic_deformation) function with `follow=:assembly` and other appropriate arguments.

````@example cHALECheckedRollManeuver
@suppress begin #hide
plot_dynamic_deformation(dynamicProblem,follow=:assembly,view=(30,15),plotDistLoads=false,plotFrequency=4,fps=24,plotLimits=([-20,20],[-20,20],[-20,20]),save=true,savePath="/docs/build/literate/cHALECheckedRollManeuver_motion.gif")
nothing #hide
end #hide
````

![](cHALECheckedRollManeuver_motion.gif)

````@example cHALECheckedRollManeuver
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

