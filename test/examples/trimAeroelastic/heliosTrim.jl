# # Trimming a flying-wing HALE
# This example illustrates how to set up a trim analysis of an aircraft in free flight. For that we take a model of the [Helios aircraft](https://en.wikipedia.org/wiki/AeroVironment_Helios_Prototype), a 75 meters wingspan flying-wing HALE. This aircraft model was described in [Patil and Hodges](https://doi.org/10.2514/1.17640):

#md # ![](../assets/helios.png)
#md # *Flying-wing HALE model geometry* by [Patil and Hodges](https://doi.org/10.2514/1.17640)

#md # ![](../assets/heliosProps.png)
#md # *Flying-wing properties* by [Patil and Hodges](https://doi.org/10.2514/1.17640)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/trimAeroelastic/heliosTrim.jl).

# ### Problem setup
# Let's begin by setting up the variables of our problem. We select the airfoil for the wing, and the optional flags to include the vertical beam pods and of setting the payload on the wing. Our set up matches that of [Patil and Hodges](https://doi.org/10.2514/1.17640). We will investigate the longitudinal trim of the aircraft at both rigid and elastic (nominal) configurations, using the payload at the center station as the varying parameter.
using AeroBeams

## Flag to include beam pods 
beamPods = true

## Flag to set payload on wing
payloadOnWing = false

## Aerodynamic solver
aeroSolver = Indicial()

## System solver 
NR = create_NewtonRaphson(absTol=1e-6,relTol=1e-6,trimSolverMethod=:LM)

## Airspeed [m/s]
U = feet2meters(40)

## Stiffness factor range
λRange = [1,50]

## Payload range [lb]
PRange = collect(0:20:500)

## Initialize outputs
problem = Array{TrimProblem}(undef,length(λRange),length(PRange))
trimAoA = Array{Float64}(undef,length(λRange),length(PRange))
trimThrust = Array{Float64}(undef,length(λRange),length(PRange))
trimδ = Array{Float64}(undef,length(λRange),length(PRange))
#md nothing #hide

# ### Solving the problem
# At this point we can sweep the stiffness factor and payload vectors to find the trim solution at each value.
## Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    ## Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress #src
        println("Trimming for λ = $λ payload = $P lb") #src
        ## Model
        helios,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,payloadOnWing=payloadOnWing,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        ## Initial guess for the trim problem
        x0 = j==1 ? zeros(0) : problem[i,j-1].systemArrays.x
        ## Create and solve trim problem
        problem[i,j] = create_TrimProblem(model=helios,systemSolver=NR,x0=x0)
        solve!(problem[i,j])
        ## Trim results
        trimSol = problem[i,j].trimSolution
        aeroVars = problem[i,j].loadstepHistory.aeroVariables
        trimAoA[i,j] = rad2deg(aeroVars[end][midSpanElem].unsteadyKin.αₑ)
        trimThrust[i,j] = trimSol.trimLoads[1]
        trimδ[i,j] = rad2deg(trimSol.trimδ[1])
        println("AoA = $(round(trimAoA[i,j], digits=2)), T = $(round(trimThrust[i,j], digits=2)), δ = $(round(trimδ[i,j], digits=2))") #src
    end
end
#md nothing #hide

# ### Post-processing
# The first post-processing step is to load the reference data by [Patil and Hodges](https://doi.org/10.2514/1.17640).
## Load reference data
αFlexibleRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_AoA_flexible.txt")
αRigidRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_AoA_rigid.txt")
δFlexibleRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_delta_flexible.txt")
δRigidRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_delta_rigid.txt")
TRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_thrust.txt")
#md nothing #hide

#md # We now plot the trim outputs as functions of the payload, for both rigid and elastic configurations. The agreement with the reference is very good. Notice that the thrust (per motor) remains almost constant, for both rigid and flexible aircraft, because the drag slope with flap deflection is assumed null (``c_{d_\delta} = 0``) and induced drag is neglected, so that only the parasite drag (``c_{d_0}``) component is present.
## Plot configurations
#md ENV["GKSwstype"] = "100" #hide
#md colors = cgrad(:rainbow, length(λRange), categorical=true)
#md labels = ["Flexible" "Rigid"]
#md lw = 2
#md ms = 3
#md msw = 0
#md ts = 10
#md fs = 16
#md gr()
#md nothing #hide

## Trim root angle of attack
#md plt1 = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,5], tickfont=font(ts), guidefont=font(fs), legend=:bottomleft, legendfontsize=12)
#md plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
#md for (i,λ) in enumerate(λRange)
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (i,λ) in enumerate(λRange)
#md     plot!(PRange, trimAoA[i,:], c=colors[i], lw=lw, label=false)
#md     if i==1
#md         scatter!(αFlexibleRef[1,:], αFlexibleRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
#md     else
#md         scatter!(αRigidRef[1,:], αRigidRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
#md     end
#md end
#md savefig("heliosTrim_AoA.svg") #hide

## Trim elevator deflection
#md plt2 = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,8], tickfont=font(ts), guidefont=font(fs))
#md for (i,λ) in enumerate(λRange)
#md     plot!(PRange, trimδ[i,:], c=colors[i], lw=lw, label=false)
#md     if i==1
#md         scatter!(δFlexibleRef[1,:], δFlexibleRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
#md     else
#md         scatter!(δRigidRef[1,:], δRigidRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
#md     end
#md end
#md savefig("heliosTrim_delta.svg") #hide

## Trim thrust per motor
#md plt3 = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500], ylims=[0,40], tickfont=font(ts), guidefont=font(fs))
#md for (i,λ) in enumerate(λRange)
#md     plot!(PRange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
#md     scatter!(TRef[1,:], TRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
#md end
#md savefig("heliosTrim_thrust.svg") #hide
#md nothing #hide

#md # ![](heliosTrim_AoA.svg)
#md # ![](heliosTrim_delta.svg)
#md # ![](heliosTrim_thrust.svg)

## Plots
# Set paths
relPath = "/test/outputs/figures/heliosTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape of the flexible aircraft at highest payload
deformationPlot = plot_steady_deformation(problem[1,end],view=(30,30),save=true,savePath=string(relPath,"/heliosTrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = cgrad(:rainbow, length(λRange), categorical=true)
labels = ["Flexible" "Rigid"]
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
gr()

# Trim root angle of attack
plt1 = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,5], tickfont=font(ts), guidefont=font(fs), legend=:bottomleft, legendfontsize=12)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
for (i,λ) in enumerate(λRange)
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimAoA[i,:], c=colors[i], lw=lw, label=false)
    if i==1
        scatter!(αFlexibleRef[1,:], αFlexibleRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    else
        scatter!(αRigidRef[1,:], αRigidRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/heliosTrim_AoA.pdf"))

# Trim thrust force
plt2 = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500], ylims=[0,40], tickfont=font(ts), guidefont=font(fs))
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
    scatter!(TRef[1,:], TRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosTrim_thrust.pdf"))

# Trim elevator deflection
plt3 = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,8], tickfont=font(ts), guidefont=font(fs))
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimδ[i,:], c=colors[i], lw=lw, label=false)
    if i==1
        scatter!(δFlexibleRef[1,:], δFlexibleRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    else
        scatter!(δRigidRef[1,:], δRigidRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/heliosTrim_delta.pdf"))

println("Finished heliosTrim.jl") #src