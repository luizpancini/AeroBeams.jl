using AeroBeams

# Flight altitude [m]
h = feet2meters(10e3)

# Airspeed [m/s]
U = 75

# Aerodynamic solver
aeroSolver = Indicial()

# Set system solver for trim problem
NRtrim = create_NewtonRaphson(trimSolverMethod=:pinv)

# Model for trim problem
BWBtrim,_ = create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=BWBtrim,systemSolver=NRtrim)
solve!(trimProblem)

# Extract trim variables
trimThrust = trimProblem.trimSolution.trimLoads[1]
trimδ = trimProblem.trimSolution.trimδ[1]
println("Trim variables: T = $(round(trimThrust, digits=2)) N, δ = $(round(rad2deg(trimδ), digits=2)) deg")

# Set checked elevator deflection profile
Δδ = deg2rad(-2)
tδinit = 1
tδramp = 1
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
BWBdynamic,_ = create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElev=δ,thrust=trimThrust,updateAirfoilParameters=false)

# Time variables
Δt = 1e-2
tf = 60

# System solver for dynamic problem
NRdyn = create_NLsolver(method=:newton,linesearch=false)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=BWBdynamic,x0=trimProblem.systemArrays.x[BWBtrim.assemblyData.rangeNominal],finalTime=tf,Δt=Δt,systemSolver=NRdyn)
@time solve!(dynamicProblem)

# Unpack numerical solution
sol = dynamicProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
rootAoA = rad2deg.([(aeroVars[i][11].unsteadyKin.αₑ+aeroVars[i][12].unsteadyKin.αₑ)/2 for i in 1:length(t)])
Δu1 = [states[i][11].nodal.u2[1] for i in 1:length(t)] .- states[1][11].nodal.u2[1]
Δu2 = [states[i][11].nodal.u2[2] for i in 1:length(t)] .- states[1][11].nodal.u2[2]
Δu3 = [states[i][11].nodal.u2[3] for i in 1:length(t)] .- states[1][11].nodal.u2[3]

## Plots
# Set paths
relPath = "/test/outputs/figures/BWBcheckedPitchManeuver"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:assembly,view=(30,15),plotDistLoads=false,plotFrequency=10,plotLimits=([-5,5],[-5,5],[-5,5]),save=true,savePath=string(relPath,"/BWBcheckedPitchManeuver_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Altitude
plt1 = plot(xlabel="Time [s]", ylabel="Altitude [m]")
plot!(t, Δu3, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBcheckedPitchManeuver_altitude.pdf"))

# Root AoA
plt2 = plot(xlabel="Time [s]", ylabel="Root angle of attack [deg]")
plot!(t, rootAoA, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/BWBcheckedPitchManeuver_rootAoA.pdf"))

println("Finished BWBcheckedPitchManeuver.jl")