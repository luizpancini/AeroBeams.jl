using AeroBeams, DelimitedFiles

# Flare angle range [rad]
ΛRange = π/180*[10,20,30]

# Root pitch angle range
θRange = π/180*vcat(-20:1:30)

# Airspeed
U = 22

# Solution method for hinge constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = false

# Spring stiffness
kSpring = 1e-4
kIPBendingHinge = 1e-1

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution, especially at lower airspeeds)
hasTipCorrection = true
tipLossDecayFactor = 8

# Initialize model
HealySideslipFFWTsteadyFlareRangeAoARangeCoast = create_HealySideslipFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,flareAngle=0,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,pitchAngle=0,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# System solver
σ0 = 1
maxIter = 200
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Fixed properties of the model
elemNodes = vcat([vcat(HealySideslipFFWTsteadyFlareRangeAoARangeCoast.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [HealySideslipFFWTsteadyFlareRangeAoARangeCoast.r_n[n][1] for n in elemNodes]
r_n2 = [HealySideslipFFWTsteadyFlareRangeAoARangeCoast.r_n[n][2] for n in elemNodes]
r_n3 = [HealySideslipFFWTsteadyFlareRangeAoARangeCoast.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(HealySideslipFFWTsteadyFlareRangeAoARangeCoast.beams[1].elements[e].x1_n1,HealySideslipFFWTsteadyFlareRangeAoARangeCoast.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [HealySideslipFFWTsteadyFlareRangeAoARangeCoast.beams[1].elements[e].x1 for e in 1:15]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
u2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
u3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
p1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
p2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
p3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
M2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
α = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
cn = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
pHinge = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
ϕHinge = Array{Float64}(undef,length(ΛRange),length(θRange))
problem = Array{SteadyProblem}(undef,length(ΛRange),length(θRange))

# Loop flare angle
for (i,Λ) in enumerate(ΛRange)
    # Loop root pitch angle
    for (j,θ) in enumerate(θRange)
        println("Solving for Λ = $(round(Int,Λ*180/π)) deg, θ = $(round(Int,θ*180/π)) deg")
        # Update model
        model = create_HealySideslipFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,flareAngle=Λ,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
        # Set initial guess solution as the one from previous pitch angle
        x0 = (j==1) ? zeros(0) : problem[i,j-1].x
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,systemSolver=NR,x0=x0)
        solve!(problem[i,j])
        # Get outputs
        u1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[1],problem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...)
        u2[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[2],problem[i,j].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:15]...)
        u3[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[3],problem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
        p1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].p_n1[1],problem[i,j].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
        p2[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].p_n1[2],problem[i,j].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
        p3[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].p_n1[3],problem[i,j].nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...)
        M2[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].M_n1[2],problem[i,j].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)
        α[i,j] = [problem[i,j].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
        cn[i,j] = [problem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]
        pHinge[i,j] = problem[i,j].model.hingeAxisConstraints[1].pH
        ϕHinge[i,j] = problem[i,j].model.hingeAxisConstraints[1].ϕ*180/π
    end
end

# Load reference data
flare10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare10_exp.txt")
flare10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare10_num.txt")
flare20_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare20_exp.txt")
flare20_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare20_num.txt")
flare30_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare30_exp.txt")
flare30_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare30_num.txt")

println("Finished HealySideslipFFWTsteadyFlareRangeAoARangeCoast.jl")