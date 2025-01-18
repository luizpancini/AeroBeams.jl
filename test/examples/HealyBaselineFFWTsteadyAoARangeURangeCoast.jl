using AeroBeams, DelimitedFiles

# Hinge configuration
hingeConfiguration = "free"

# Root pitch angle range
θRange = π/180*[2.5; 5.0; 7.5]

# Airspeed Range
URange = vcat(0:1:40)

# Flare angle [rad]
Λ = 15*π/180

# Gravity
g = 9.80665

# Discretization
nElementsInner = 16
nElementsFFWT = 4

# Tip loss options (assumed, since Healy's analysis uses DLM for aerodynamic)
withTipCorrection = true
tipLossDecayFactor = 12

# Initialize model
HealyBaselineFFWTsteadyAoARangeURangeCoast = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,flareAngle=Λ,pitchAngle=0,withTipCorrection=withTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# System solver
σ0 = 1
maxIter = 200
relTol = 1e-6
ΔλRelaxFactor = 1/2
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Fixed properties of the model
hingeNode = nElementsInner+1
L = HealyBaselineFFWTsteadyAoARangeURangeCoast.beams[1].length
chords = [HealyBaselineFFWTsteadyAoARangeURangeCoast.beams[1].elements[e].aero.c for e in 1:15]
Δℓs = [HealyBaselineFFWTsteadyAoARangeURangeCoast.beams[1].elements[e].Δℓ for e in 1:15]
ρ = HealyBaselineFFWTsteadyAoARangeURangeCoast.atmosphere.ρ

# Initialize outputs
lift = Array{Float64}(undef,length(θRange),length(URange))
u3Hinge = Array{Float64}(undef,length(θRange),length(URange))
M2root = Array{Float64}(undef,length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(θRange),length(URange))
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))

# Loop root pitch angle
for (i,θ) in enumerate(θRange)
    # Loop airspeed
    for (j,U) in enumerate(URange)
        println("Solving for θ = $(round(θ*180/π,digits=1)) deg, U = $U m/s")
        # Update model
        model = create_HealyBaselineFFWT(hingeConfiguration=hingeConfiguration,flareAngle=Λ,airspeed=U,pitchAngle=θ,withTipCorrection=withTipCorrection,tipLossDecayFactor=tipLossDecayFactor,g=g,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
        # Set initial guess solution as the one from previous airspeed
        x0 = (j==1) ? zeros(0) : problem[i,j-1].x
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,systemSolver=NR,x0=x0)
        solve!(problem[i,j])
        # Get outputs
        u3Hinge[i,j] = problem[i,j].nodalStatesOverσ[end][nElementsInner].u_n2[3]
        M2root[i,j] = problem[i,j].nodalStatesOverσ[end][1].M_n1[2]
        α = [problem[i,j].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
        cn = [problem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]
        ct = [problem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.ct for e in 1:15]
        cl = @. cn*cos(α) + ct*sin(α)
        lift[i,j] = sum(1/2*ρ*U^2*chords.*Δℓs.*cl)
        ϕHinge[i,j] = problem[i,j].model.hingeAxisConstraints[1].ϕ*180/π
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

println("Finished HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")