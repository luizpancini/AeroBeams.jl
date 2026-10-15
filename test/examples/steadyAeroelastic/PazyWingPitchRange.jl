# # Steady aeroelastic analysis of the Pazy wing
# This example illustrates how to set up a steady aeroelastic analysis, using the Technion's Pazy wing benchmark. The sectional properties of the wing's spar and the aerodynamic tip loss function were defined by [Riso and Cesnik](https://doi.org/10.2514/6.2022-2313). The data is publicly available at [https://github.com/UM-A2SRL/AePW3-LDWG](https://github.com/UM-A2SRL/AePW3-LDWG).

#md # ![](../assets/pazyWT.png)

#md # *Pazy wing in the wind tunnel* by [Avin et al.](https://doi.org/10.2514/1.J060621)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/steadyAeroelastic/PazyWingPitchRange.jl).

# ### Problem setup
# Let's begin by setting the variables of our problem. In this example we will analyze the displacements and twist of the clamped wing under several combinations of root pitch angle and airspeed, which are defined by the arrays `θRange` and `URange`. The wing was tested with a balance mass of 0.01 kg attached to its tip store, at a position of 40 mm behind the trailing-edge. Notice that we bring into scope some fixed geometrical and discretization properties of the wing's beam through the function [`geometrical_properties_Pazy`](@ref geometrical_properties_Pazy).
using AeroBeams

## Aerodynamic solver
aeroSolver = Indicial()

## Airfoil section
airfoil = deepcopy(FlatPlate)

## Flag for upright position
upright = true

## Gravity
g = 0

## Flag for small angles approximation
smallAngles = true

## Fixed geometrical and discretization properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy()

## Tip mass (for test 1, 0.01 kg, 40 mm behind the trailing-edge)
tipMass = 0.01
ηtipMass = [0; -chord*(1-normSparPos)-0.04; 0]

## Root pitch angle [rad] and airspeed [m/s] ranges
θRange = deg2rad.([3,5,7])
URange = vcat(1e-3,1:1:100)

## Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
tip_IP = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))
tip_AoA = Array{Float64}(undef,length(θRange),length(URange))
#md nothing #hide

# ### Solving the problem
# In the following loops, we create new model instances with the combination of pitch angle and airspeed, create and solve the steady problem, and then extract the outputs of interest. The model creation process is streamlined with the function [`create_Pazy`](@ref create_Pazy), taking the appropriate inputs. Notice that the twist angle is computed as the arcsine of the difference in "height" between the leading edge and trailing edge of the wingtip (throught the rotation of the local basis vector pointing in the x2-direction).
## Sweep root angle
for (i,θ) in enumerate(θRange)
    ## Initialize guess solution
    x0 = zeros(0)
    ## Sweep airspeed
    for (j,U) in enumerate(URange)
        ## Display progress #src
        println("Solving for θ = $(round(rad2deg(θ))) deg, U = $U m/s") #src
        ## Update model
        PazyWingPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,pitchAngle=θ,airspeed=U,gravityVector=[0,0,-g],tipMass=tipMass,ηtipMass=ηtipMass,smallAngles=smallAngles)
        ## Create and solve problem
        problem[i,j] = create_SteadyProblem(model=PazyWingPitchRange,x0=x0)
        solve!(problem[i,j])
        ## Get tip twist, AoA, IP and OOP displacement at beam reference line
        tipElemFinalAero = problem[i,j].loadstepHistory.aeroVariables[end][nElem]
        tipNodeFinalState = problem[i,j].loadstepHistory.states[end][nElem].nodal
        tip_twist[i,j] = twist_angle(tipNodeFinalState.p2_b)
        tip_OOP[i,j] = tipNodeFinalState.u2_b[3]
        tip_IP[i,j] = -tipNodeFinalState.u2_b[2]
        tip_AoA[i,j] = rad2deg(tipElemFinalAero.unsteadyKin.αₑ)
        ## Set initial guess solution for next iteration
        x0 = problem[i,j].systemArrays.x
    end
end
#md nothing #hide

# ### Post-processing
# The post-processing begins by loading the reference data.
## Load reference data
tip_u3VsU_rootPitch5_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch5_Exp.txt")
tip_u3VsU_rootPitch5_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch5_UMNAST.txt")
tip_u3VsU_rootPitch7_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch7_Exp.txt")
tip_u3VsU_rootPitch7_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch7_UMNAST.txt")
tip_thetaVsU_rootPitch5_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch5_Exp.txt")
tip_thetaVsU_rootPitch5_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch5_UMNAST.txt")
tip_thetaVsU_rootPitch7_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch7_Exp.txt")
tip_thetaVsU_rootPitch7_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch7_UMNAST.txt")
tip_u3Vsq_rootPitch3_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch3_ExpTest1.txt")
tip_u3Vsq_rootPitch5_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch5_ExpTest1.txt")
tip_u3Vsq_rootPitch7_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch7_ExpTest1.txt")
ρTest2 = 2*1050/43^2 # Estimated air density for test 2, from Table 13 of Avin et al.'s paper
#md nothing #hide

#md # We can now plot the outputs as a function of airspeed for each of the root pitch angles. The following "experimental" results were taken from Figure 33 of the paper by [Avin et al.](https://doi.org/10.2514/1.J060621). Note that the experimental twist angle was actually estimated by [Riso and Cesnik](https://arc.aiaa.org/doi/abs/10.2514/6.2023-0759) using the diffence between the leading- and trailing-edge out-of-plane displacements of [Avin et al.](https://doi.org/10.2514/1.J060621). The correlation with the experimental data and the reference numerical solution is very good.
#md using Suppressor #hide
#md colors = cgrad(:rainbow,  length(θRange), categorical=true)
#md gr()
#md ENV["GKSwstype"] = "100" #hide

## Tip OOP displacement vs. airspeed
#md plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP displacement [% semispan]", xlims=[0,60], ylims=[0,50], legend=:topleft)
#md plot!([NaN], [NaN], c=:black, lw=2, ls=:solid, label="AeroBeams")
#md plot!([NaN], [NaN], c=:black, lw=2, ls=:dashdot, label="UM/NAST")
#md scatter!([NaN], [NaN], c=:black, ms=4, label="Exp. (test 1)")
#md scatter!([NaN], [NaN], c=:black, shape=:diamond, ms=4, label="Exp. (test 2)")
#md for (i,θ) in enumerate(θRange)
#md     plot!(URange, tip_OOP[i,:]/L*100, c=colors[i], lw=2, ls=:solid, label=string("\$\\alpha_r=",round(Int,rad2deg(θ)),"^\\circ\$"))
#md     if i==1
#md         scatter!(sqrt.(2*tip_u3Vsq_rootPitch3_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch3_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=4, msw=0, label=false)
#md     elseif i==2
#md         plot!(tip_u3VsU_rootPitch5_UMNAST[1,:], tip_u3VsU_rootPitch5_UMNAST[2,:], lw=2, ls=:dashdot, c=colors[i], label=false)
#md         scatter!(tip_u3VsU_rootPitch5_Exp[1,:], tip_u3VsU_rootPitch5_Exp[2,:], mc=colors[i], ms=4, msw=0, label=false)
#md         scatter!(sqrt.(2*tip_u3Vsq_rootPitch5_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch5_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=4, msw=0, label=false)
#md     elseif i==3
#md         plot!(tip_u3VsU_rootPitch7_UMNAST[1,:], tip_u3VsU_rootPitch7_UMNAST[2,:], lw=2, ls=:dashdot, c=colors[i], label=false)
#md         scatter!(tip_u3VsU_rootPitch7_Exp[1,:], tip_u3VsU_rootPitch7_Exp[2,:], mc=colors[i], ms=4, msw=0, label=false)
#md         scatter!(sqrt.(2*tip_u3Vsq_rootPitch7_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch7_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=4, msw=0, label=false)
#md     end
#md end
#md savefig("PazyWingPitchRange_tipOOP.svg") #hide
#md nothing #hide

## Tip twist vs. airspeed
#md plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,60], legend=:topleft)
#md for (i,θ) in enumerate(θRange)
#md     plot!(URange, tip_twist[i,:], c=colors[i], lw=2, label=false)
#md     if i==2
#md         plot!(tip_thetaVsU_rootPitch5_UMNAST[1,:], tip_thetaVsU_rootPitch5_UMNAST[2,:], lw=2, ls=:dashdot, c=colors[i], label=false)
#md         scatter!(tip_thetaVsU_rootPitch5_Exp[1,:], tip_thetaVsU_rootPitch5_Exp[2,:], mc=colors[i], ms=4, msw=0, label=false)
#md     elseif i==3
#md         plot!(tip_thetaVsU_rootPitch7_UMNAST[1,:], tip_thetaVsU_rootPitch7_UMNAST[2,:], lw=2, ls=:dashdot, c=colors[i], label=false)
#md         scatter!(tip_thetaVsU_rootPitch7_Exp[1,:], tip_thetaVsU_rootPitch7_Exp[2,:], mc=colors[i], ms=4, msw=0, label=false)
#md     end
#md end
#md savefig("PazyWingPitchRange_tipTwist.svg") #hide
#md nothing #hide

## Tip AoA vs. airspeed
#md plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip AoA [deg]", xlims=[0,60], legend=:bottomright)
#md for (i,θ) in enumerate(θRange)
#md     plot!(URange, tip_AoA[i,:], c=colors[i], lw=2, label=false)
#md end
#md savefig("PazyWingPitchRange_tipAoA.svg") #hide
#md nothing #hide

## Tip in-plane displacement vs. airspeed 
#md plt4 = plot(xlabel="Airspeed [m/s]", ylabel="Tip IP displacement [% semispan]", xlims=[0,60], legend=:topleft)
#md for (i,θ) in enumerate(θRange)
#md     plot!(URange, tip_IP[i,:]/L*100, c=colors[i], lw=2, label=false)
#md end
#md savefig("PazyWingPitchRange_tipIP.svg") #hide
#md nothing #hide

#md # ![](PazyWingPitchRange_tipOOP.svg)
#md # ![](PazyWingPitchRange_tipTwist.svg)
#md # ![](PazyWingPitchRange_tipAoA.svg)
#md # ![](PazyWingPitchRange_tipIP.svg)

# Let's leverage the built-in function [`plot_steady_outputs`](@ref plot_steady_outputs) with the appropriate inputs to plot the lift coefficient over the span of the wing.
## Lift coefficient over span for the last problem (θ = θRange[end] and U = URange[end])
#md @suppress begin #hide
#md plot_steady_outputs(problem[end],outputs=["cl"],colorScheme=:grays,lw=2)
#md savefig("PazyWingPitchRange_cl.svg") #hide
#md nothing #hide
#md end #hide
#md # ![](PazyWingPitchRange_cl.svg)
#md nothing #hide

#md # Finally, we can visualize the deformed state of the wing using the [`plot_steady_deformation`](@ref plot_steady_deformation) with the appropriate arguments. Notice that by default, gravitational loads are plotted in yellow, and aerodynamic loads in green.
## Deformed state of last problem (θ = θRange[end] and U = URange[end])
#md @suppress begin #hide
#md deformationPlot = plot_steady_deformation(problem[end],view=(45,30),plotBCs=true,plotDistLoads=true)
#md savefig("PazyWingPitchRange_deformation.svg") #hide
#md end #hide
#md nothing #hide
#md # ![](PazyWingPitchRange_deformation.svg)
#md nothing #hide

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed state of last problem
deformationPlot = plot_steady_deformation(problem[end],view=(45,30),save=true,savePath="/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_deformation.pdf")
display(deformationPlot)

# Plot configurations
colors = cgrad(:rainbow,length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
msw = 0
gr()

# Tip OOP displacement vs. airspeed for root several pitch angles 
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,60], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dashdot, label="UM/NAST (Exp. loss) - Riso & Cesnik (2022,2023)")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Exp. (test 1) - Avin et al. (2022)")
scatter!([NaN], [NaN], c=:black, shape=:diamond, ms=ms, label="Exp. (test 2) - Avin et al. (2022)")
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_OOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label=string("\$\\alpha_r=",round(Int,rad2deg(θ)),"^\\circ\$"))
    if i==1
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch3_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch3_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif i==2
        plot!(tip_u3VsU_rootPitch5_UMNAST[1,:], tip_u3VsU_rootPitch5_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_u3VsU_rootPitch5_Exp[1,:], tip_u3VsU_rootPitch5_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch5_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch5_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif i==3
        plot!(tip_u3VsU_rootPitch7_UMNAST[1,:], tip_u3VsU_rootPitch7_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_u3VsU_rootPitch7_Exp[1,:], tip_u3VsU_rootPitch7_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch7_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch7_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/PazyWingPitchRange_tipOOP.pdf"))

# Tip twist vs. airspeed for root several pitch angles 
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_twist[i,:], c=colors[i], lw=lw, label=false)
    if i==2
        plot!(tip_thetaVsU_rootPitch5_UMNAST[1,:], tip_thetaVsU_rootPitch5_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_thetaVsU_rootPitch5_Exp[1,:], tip_thetaVsU_rootPitch5_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
    elseif i==3
        plot!(tip_thetaVsU_rootPitch7_UMNAST[1,:], tip_thetaVsU_rootPitch7_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_thetaVsU_rootPitch7_Exp[1,:], tip_thetaVsU_rootPitch7_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/PazyWingPitchRange_tipTwist.pdf"))

# Tip AoA vs. airspeed for root several pitch angles 
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip AoA [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_AoA[i,:], c=colors[i], lw=lw, label=false)
end
display(plt3)
savefig(string(absPath,"/PazyWingPitchRange_tipAoA.pdf"))

# Tip in-plane displacement vs. airspeed for root several pitch angles 
plt4 = plot(xlabel="Airspeed [m/s]", ylabel="Tip IP disp. [% semispan]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_IP[i,:]/L*100, c=colors[i], lw=lw, label=false)
end
display(plt4)
savefig(string(absPath,"/PazyWingPitchRange_tipIP.pdf"))

# Lift coefficient over span for the last solved problem
plot_steady_outputs(problem[end],outputs=["cl","κ2","F3"],colorScheme=:grays,lw=lw,save=true,saveFolder=string(relPath,"/"))

println("Finished PazyWingPitchRange.jl") #src