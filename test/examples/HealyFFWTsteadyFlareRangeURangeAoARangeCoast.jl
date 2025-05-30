using AeroBeams, DelimitedFiles

# Flare angle range [rad]
ΛRange = π/180*[10,20,30]

# Root pitch angle range
θRange = π/180*[0,5,10]

# Airspeed range
URange = collect(10:1:35)

# Spring stiffness
kSpring = 1e-4

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution, especially at lower airspeeds)
hasTipCorrection = true
tipLossDecayFactor = 10

# Initialize model
HealyFFWTsteadyFlareRangeURangeAoARangeCoast = create_HealyFFWT(flareAngle=0,kSpring=kSpring,pitchAngle=0,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# System solver
σ0 = 1
maxIter = 200
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Fixed properties of the model
elemNodes = vcat([vcat(HealyFFWTsteadyFlareRangeURangeAoARangeCoast.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [HealyFFWTsteadyFlareRangeURangeAoARangeCoast.r_n[n][1] for n in elemNodes]
r_n2 = [HealyFFWTsteadyFlareRangeURangeAoARangeCoast.r_n[n][2] for n in elemNodes]
r_n3 = [HealyFFWTsteadyFlareRangeURangeAoARangeCoast.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(HealyFFWTsteadyFlareRangeURangeAoARangeCoast.beams[1].elements[e].x1_n1,HealyFFWTsteadyFlareRangeURangeAoARangeCoast.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [HealyFFWTsteadyFlareRangeURangeAoARangeCoast.beams[1].elements[e].x1 for e in 1:15]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
u2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
u3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
p1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
p2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
p3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
M2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
α = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
pHinge = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(ΛRange),length(θRange),length(URange))
problem = Array{SteadyProblem}(undef,length(ΛRange),length(θRange),length(URange))

# Loop flare angle
for (i,Λ) in enumerate(ΛRange)
    # Loop root pitch angle
    for (j,θ) in enumerate(θRange)
        # Loop airspeed
        for (k,U) in enumerate(URange)
            println("Solving for Λ = $(round(Int,Λ*180/π)) deg, θ = $(round(Int,θ*180/π)) deg, U = $U m/s")
            # Update model
            model = create_HealyFFWT(flareAngle=Λ,kSpring=kSpring,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
            # Set initial guess solution as the one from previous airspeed
            x0 = (k==1) ? zeros(0) : problem[i,j,k-1].x
            # Create and solve problem
            problem[i,j,k] = create_SteadyProblem(model=model,systemSolver=NR,x0=x0)
            solve!(problem[i,j,k])
            # Get outputs
            u1[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].u_n1[1],problem[i,j,k].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...)
            u2[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].u_n1[2],problem[i,j,k].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:15]...)
            u3[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].u_n1[3],problem[i,j,k].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
            p1[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].p_n1[1],problem[i,j,k].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
            p2[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].p_n1[2],problem[i,j,k].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
            p3[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].p_n1[3],problem[i,j,k].nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...)
            M2[i,j,k] = vcat([vcat(problem[i,j,k].nodalStatesOverσ[end][e].M_n1[2],problem[i,j,k].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)
            α[i,j,k] = [problem[i,j,k].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
            cn[i,j,k] = [problem[i,j,k].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]
            pHinge[i,j,k] = problem[i,j,k].model.hingeAxisConstraints[1].pH
            ϕHinge[i,j,k] = problem[i,j,k].model.hingeAxisConstraints[1].ϕ*180/π
        end
    end
end

# Load reference data
flare10_aoa0_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa0_exp.txt")
flare10_aoa0_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa0_num.txt")
flare10_aoa5_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa5_exp.txt")
flare10_aoa5_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa5_num.txt")
flare10_aoa10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa10_exp.txt")
flare10_aoa10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa10_num.txt")
flare20_aoa0_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa0_exp.txt")
flare20_aoa0_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa0_num.txt")
flare20_aoa5_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa5_exp.txt")
flare20_aoa5_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa5_num.txt")
flare20_aoa10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa10_exp.txt")
flare20_aoa10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa10_num.txt")
flare30_aoa0_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa0_exp.txt")
flare30_aoa0_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa0_num.txt")
flare30_aoa5_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa5_exp.txt")
flare30_aoa5_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa5_num.txt")
flare30_aoa10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa10_exp.txt")
flare30_aoa10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa10_num.txt")

println("Finished HealyFFWTsteadyFlareRangeURangeAoARangeCoast.jl")