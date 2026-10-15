using AeroBeams

# Hinge configuration
hingeConfiguration = :free

# Root pitch angle range [rad]
θRange = deg2rad.([2.5, 5.0, 7.5])

# Airspeed range [m/s]
URange = vcat(1:1:40)

# Aerodynamic solver
aeroSolver = QuasiSteady()

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options (assumed, since Healy's analysis uses DLM for aerodynamic)
hasTipCorrection = true
tipLossDecayFactor = 20

# Initialize model
HealyBaselineFFWTsteadyAoARangeURangeCoast = create_HealyBaselineFFWT(aeroSolver=aeroSolver,hingeConfiguration=hingeConfiguration,pitchAngle=0,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# System solver
# systemSolver = create_NewtonRaphson()
systemSolver = create_NLsolver(displayStatus=false,method=:newton,linesearch=true,jcalc=:analytical,checkJac=false)

# Fixed properties of the model
hingeNode = nElementsInner+1
L = HealyBaselineFFWTsteadyAoARangeURangeCoast.beams[1].geometry.ℓ
chords = [HealyBaselineFFWTsteadyAoARangeURangeCoast.beams[1].elements[e].aero.geometry.c for e in 1:nElem]
Δℓs = [HealyBaselineFFWTsteadyAoARangeURangeCoast.beams[1].elements[e].geometry.Δℓ for e in 1:nElem]
ρ = HealyBaselineFFWTsteadyAoARangeURangeCoast.atmosphere.ρ

# Initialize outputs
lift = Array{Float64}(undef,length(θRange),length(URange))
u3Hinge = Array{Float64}(undef,length(θRange),length(URange))
M2root = Array{Float64}(undef,length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(θRange),length(URange))

# Loop root pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Loop airspeed
    for (j,U) in enumerate(URange)
        println("Solving for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Update model
        model = create_HealyBaselineFFWT(aeroSolver=aeroSolver,hingeConfiguration=hingeConfiguration,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
        # Create and solve problem
        problem = create_SteadyProblem(model=model,systemSolver=systemSolver,x0=x0)
        solve!(problem)
        # Get outputs
        state = problem.loadstepHistory.states[end]
        aeroVars = problem.loadstepHistory.aeroVariables[end]
        u3Hinge[i,j] = state[nElementsInner].nodal.u2[3]
        M2root[i,j] = state[1].nodal.M1[2]
        α = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
        cn = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
        ct = [aeroVars[e].aeroCoefficients.ct for e in 1:nElem]
        cl = @. cn*cosd(α) + ct*sind(α)
        lift[i,j] = sum(1/2*ρ*U^2*chords.*Δℓs.*cl)
        ϕHinge[i,j] = -rad2deg(problem.model.hingeAxisConstraints[1].state.ϕ)
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
end

# Load reference data
fold_aoa_25 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/fold_aoa_25.txt")
fold_aoa_50 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/fold_aoa_50.txt")
fold_aoa_75 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/fold_aoa_75.txt")

lift_aoa_25 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/lift_aoa_25.txt")
lift_aoa_50 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/lift_aoa_50.txt")
lift_aoa_75 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/lift_aoa_75.txt")

RBM_aoa_25 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/RBM_aoa_25.txt")
RBM_aoa_50 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/RBM_aoa_50.txt")
RBM_aoa_75 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/RBM_aoa_75.txt")

disp_aoa_25 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/disp_aoa_25.txt")
disp_aoa_50 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/disp_aoa_50.txt")
disp_aoa_75 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/disp_aoa_75.txt")

theta_aoa_25 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/theta_aoa_25.txt")
theta_aoa_50 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/theta_aoa_50.txt")
theta_aoa_75 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTsteadyAoARangeURangeCoast/theta_aoa_75.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTsteadyAoARangeURangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
ms = 6
msw = 0
gr()

# Hinge OOP displacement
pltu = plot(xlabel="Airspeed [m/s]", ylabel="Hinge displacement [\$\\%\$ of semispan]", xlims=[0,40], ylims=[-10,20], yticks=-10:10:20, legendfontsize=10, legend=:topleft)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        disp = disp_aoa_25
    elseif i==2
        disp = disp_aoa_50
    else
        disp = disp_aoa_75
    end
    plot!(URange, u3Hinge[i,:]/L*100, lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(rad2deg(θ),digits=1)) \\degree\$")
    plot!(disp[1,:], disp[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltu)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_u3Hinge.pdf"))

# Root bending moment
pltM = plot(xlabel="Airspeed [m/s]", ylabel="Root bending moment [N.m]", xlims=[0,40], ylims=[-10,40], yticks=-10:10:40, legendfontsize=10, legend=:topleft)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        RBM = RBM_aoa_25
    elseif i==2
        RBM = RBM_aoa_50
    else
        RBM = RBM_aoa_75
    end
    plot!(URange, -M2root[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(rad2deg(θ),digits=1)) \\degree\$")
    plot!(RBM[1,:], RBM[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltM)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_M2root.pdf"))

# Lift
pltLift = plot(xlabel="Airspeed [m/s]", ylabel="Lift [N]", xlims=[0,40], ylims=[-10,60], yticks=-10:10:60, legendfontsize=10, legend=:topleft)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        liftHealy = lift_aoa_25
    elseif i==2
        liftHealy = lift_aoa_50
    else
        liftHealy = lift_aoa_75
    end
    plot!(URange, lift[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(rad2deg(θ),digits=1)) \\degree\$")
    plot!(liftHealy[1,:], liftHealy[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltLift)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_lift.pdf"))

# Coast angle
pltϕ = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", xlims=[0,40], ylims=[-90,90], yticks=-90:30:90, legendfontsize=10, legend=:bottomright)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        fold = fold_aoa_25
    elseif i==2
        fold = fold_aoa_50
    else
        fold = fold_aoa_75
    end
    plot!(URange, ϕHinge[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(rad2deg(θ),digits=1)) \\degree\$")
    plot!(fold[1,:], fold[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltϕ)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_coastAngle.pdf"))

println("Finished HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")