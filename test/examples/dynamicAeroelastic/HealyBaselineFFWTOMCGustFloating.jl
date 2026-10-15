# # Gust response of a wing with flared folding wingtip (FFWT)
# This example illustrates how to set up a dynamic gust response analysis for a wing featuring a flared folding wingtip (FFWT). The influence of a one-minus-cosine gust frequency on the behavior of this wing model was studied numerically by [Healy's PhD Thesis](https://research-information.bris.ac.uk/ws/portalfiles/portal/388426634/Final_Copy_2024_01_23_Healy_F_PhD_Redacted.pdf).

#md # ![](../assets/HealyBaselineFFWT.png)

#md # *Top view of baseline wing model featuring a FFWT* by [Healy et al.](https://doi.org/10.2514/6.2023-0376)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating.jl).

# ### Problem setup
# Let's begin by setting the variables of our problem. In this example we will analyze the response of the coasting (fold) angle and wing root bending moment (WRBM) to a series of one-minus-cosine gusts with different frequencies (meaning different durations).
using AeroBeams

## Gust frequency range [Hz]
ωRange = [1,3,10]

## Gust maximum vertical velocity [m/s]
Ug = 0.5

## Hinge configuration
hingeConfiguration = :free

## Root pitch angle [rad]
θ = deg2rad(7.5)

## Airspeed [m/s]
U = 18

## Discretization
nElementsInner = 16
nElementsFFWT = 8
nElem = nElementsInner + nElementsFFWT

## Tip loss options (adjusted to match experimental data)
hasTipCorrection = true
tipLossDecayFactor = 17

## Time variables
Δt = 2.5e-3
tf = 2

## System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

## Initialize outputs
dynamicProblem = Array{DynamicProblem}(undef,length(ωRange))
t = Array{Vector{Float64}}(undef,length(ωRange))
M2root = Array{Vector{Float64}}(undef,length(ωRange))
ϕ = Array{Vector{Float64}}(undef,length(ωRange))
#md nothing #hide

# ### Solving the problem
# In the following loop, we create new gust and model instances for each gust frequency, create and solve the initial steady problem to use as the initial solution to the dynamic problem, which yields the time responses of the coasting angle of the FFWT (`ϕ`) and the out-of-plane bending moment at the root (`M2root`). The model creation process is streamlined with the function [`create_HealyBaselineFFWT`](@ref create_HealyBaselineFFWT), taking the appropriate inputs.
#md using Suppressor #hide
## Loop gust frequency range
for (i,ω) in enumerate(ωRange)
    ## Display progress
    display("Solving for ω = $ω Hz") #src
    ## Gust
    gust = create_OneMinusCosineGust(duration=1/ω,verticalVelocity=Ug)
    ## Model
    model = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust) 
    ## Solve steady problem for initial conditions
    steadyProblem = create_SteadyProblem(model=model,systemSolver=systemSolver)
    solve!(steadyProblem)
    ## Create and solve dynamic problem
    dynamicProblem[i] = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x)
#md     @suppress begin #hide
            solve!(dynamicProblem[i])
#md     end #hide
    ## Unpack numerical solution
    sol = dynamicProblem[i].solutionHistory
    states = sol.states
    hacStates = sol.hacStates
    t[i] = sol.t
    M2root[i] = [states[j][1].nodal.M1[2] for j in eachindex(t[i])]
    ϕ[i] = [rad2deg(hacStates[j][1].ϕ) for j in eachindex(t[i])]
end
#md nothing #hide

# ### Post-processing
# First, we choose which of Healy's reference model we will use. Then, we load the corresponding data.
## Select reference data (choose between "baseline" or "ASMVLM")
refModel = "ASMVLM"

## Load reference data
ΔWRBM_Healy_1Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTOMCGustFloating/DWRBM_"*refModel*"_1Hz.txt")
ΔWRBM_Healy_3Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTOMCGustFloating/DWRBM_"*refModel*"_3Hz.txt")
ΔWRBM_Healy_10Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTOMCGustFloating/DWRBM_"*refModel*"_10Hz.txt")
fold_Healy_1Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTOMCGustFloating/fold_"*refModel*"_1Hz.txt")
fold_Healy_3Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTOMCGustFloating/fold_"*refModel*"_3Hz.txt")
fold_Healy_10Hz = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTOMCGustFloating/fold_"*refModel*"_10Hz.txt")

#md nothing #hide

#md # We can now plot the increments in the fold angle and root bending moment of the wing for each of the gusts. The following reference results were taken from Fig. 3.34 of [Healy's PhD Thesis](https://research-information.bris.ac.uk/ws/portalfiles/portal/388426634/Final_Copy_2024_01_23_Healy_F_PhD_Redacted.pdf). The evolution of the fold angle and root bending moment is closely matched for the first few peaks, but then the correlation with the reference deteriorates somewhat: AeroBeams predicts oscillations with less damping and with a different frequency content. Nevertheless, the maximum loads, which are arguably the most important outputs, are captured with reasonable accuracy.
#md gr()
#md ENV["GKSwstype"] = "100" #hide
#md colors = cgrad(:rainbow, length(ωRange), categorical=true)
#md lw = 2

## Fold angle increment
#md plt_ϕ = plot(xlabel="Time [s]", ylabel="Fold angle increment [deg]", ylims=[-20,20], tickfont=font(10), guidefont=font(16), legendfontsize=12, legend=:topright)
#md plot!([NaN],[NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
#md plot!([NaN],[NaN], c=:black, lw=lw, ls=:dashdot, label="Healy (2023)")
#md for (i,ω) in enumerate(ωRange)
#md     plot!(t[i], -(ϕ[i] .- ϕ[i][1]), c=colors[i], lw=lw, ls=:solid, label="\$\\omega_g = $ω\$ Hz")
#md     if ω==1
#md         plot!(fold_Healy_1Hz[1,:],fold_Healy_1Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
#md     elseif ω==3
#md         plot!(fold_Healy_3Hz[1,:],fold_Healy_3Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
#md     else ω==10
#md         plot!(fold_Healy_10Hz[1,:],fold_Healy_10Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
#md     end
#md end
#md savefig("HealyBaselineFFWTOMCGustFloating_phi.svg")

## Root OOP bending moment increment
#md plt_ΔWRBM = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", tickfont=font(10), guidefont=font(16))
#md for (i,ω) in enumerate(ωRange)
#md     plot!(t[i], -(M2root[i] .- M2root[i][1]), c=colors[i], lw=lw, ls=:solid, label=false)
#md     if ω==1
#md         plot!(ΔWRBM_Healy_1Hz[1,:],ΔWRBM_Healy_1Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
#md     elseif ω==3
#md         plot!(ΔWRBM_Healy_3Hz[1,:],ΔWRBM_Healy_3Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
#md     elseif ω==10
#md         plot!(ΔWRBM_Healy_10Hz[1,:],ΔWRBM_Healy_10Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
#md     end
#md end
#md savefig("HealyBaselineFFWTOMCGustFloating_DWRBM.svg")

#md nothing #hide

#md # ![](HealyBaselineFFWTOMCGustFloating_phi.svg)
#md # ![](HealyBaselineFFWTOMCGustFloating_DWRBM.svg)

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation for last dynamic problem
L = dynamicProblem[end].model.beams[1].geometry.ℓ
anim = plot_dynamic_deformations(dynamicProblem,follow=:basisA,plotFrequency=8,fps=30,plotDistLoads=false,view=(30,15),plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),save=true,savePath=string(relPath,"/HealyBaselineFFWTOMCGustFloating_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
colors = cgrad(:rainbow, length(ωRange), categorical=true)
lw = 2
gr()

# Fold angle increment
plt_ϕ = plot(xlabel="Time [s]", ylabel="Fold angle increment [deg]", ylims=[-20,20], tickfont=font(10), guidefont=font(16), legendfontsize=12, legend=:topright)
plot!([NaN],[NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN],[NaN], c=:black, lw=lw, ls=:dashdot, label="Healy (2023)")
for (i,ω) in enumerate(ωRange)
    plot!(t[i], -(ϕ[i] .- ϕ[i][1]), c=colors[i], lw=lw, ls=:solid, label="\$\\omega_g = $ω\$ Hz")
    if i==1
        plot!(fold_Healy_1Hz[1,:],fold_Healy_1Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif i==2
        plot!(fold_Healy_3Hz[1,:],fold_Healy_3Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    else i==3
        plot!(fold_Healy_10Hz[1,:],fold_Healy_10Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    end
end
display(plt_ϕ)
savefig(string(absPath,"/HealyBaselineFFWTOMCGustFloating_fold_angle_increment.pdf"))

plt_ΔWRBM = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", tickfont=font(10), guidefont=font(16))
for (i,ω) in enumerate(ωRange)
    plot!(t[i], -(M2root[i] .- M2root[i][1]), c=colors[i], lw=lw, ls=:solid, label=false)
    if i==1
        plot!(ΔWRBM_Healy_1Hz[1,:],ΔWRBM_Healy_1Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif i==2
        plot!(ΔWRBM_Healy_3Hz[1,:],ΔWRBM_Healy_3Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    elseif i==3
        plot!(ΔWRBM_Healy_10Hz[1,:],ΔWRBM_Healy_10Hz[2,:], c=colors[i], lw=lw, ls=:dashdot, label=false)
    end
end
display(plt_ΔWRBM)
savefig(string(absPath,"/HealyBaselineFFWTOMCGustFloating_DWRBM.pdf"))

println("Finished HealyBaselineFFWTOMCGustFloating.jl") #src