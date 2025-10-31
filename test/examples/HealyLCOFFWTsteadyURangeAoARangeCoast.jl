using AeroBeams, DelimitedFiles

# Wingtip configuration
wingtipTab = false

# Root pitch angle range
θRange = π/180*[2.5, 5, 7.5, 10]

# Airspeed range
URange = collect(5:0.5:30)

# Spring stiffness
kSpring = 1e-4
kIPBendingHinge = 1e2

# Discretization
nElementsInner = 8
nElementsOuter = 3
nElementsFFWT = 4

# Tip loss options (the value of tipLossDecayFactor is assumed to match the experimental results, since it strongly influences the solution)
hasTipCorrection = true
tipLossDecayFactor = nothing

# Initialize model
HealyLCOFFWTsteadyURangeAoARangeCoast = create_HealyLCOFFWT(wingtipTab=wingtipTab,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Fixed properties of the model
nElem = nElementsInner + nElementsOuter + nElementsFFWT
elemNodes = vcat([vcat(HealyLCOFFWTsteadyURangeAoARangeCoast.elements[e].nodesGlobalID) for e in 1:nElem]...)
r_n1 = [HealyLCOFFWTsteadyURangeAoARangeCoast.r_n[n][1] for n in elemNodes]
r_n2 = [HealyLCOFFWTsteadyURangeAoARangeCoast.r_n[n][2] for n in elemNodes]
r_n3 = [HealyLCOFFWTsteadyURangeAoARangeCoast.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(HealyLCOFFWTsteadyURangeAoARangeCoast.beams[1].elements[e].x1_n1,HealyLCOFFWTsteadyURangeAoARangeCoast.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
x1_e = [HealyLCOFFWTsteadyURangeAoARangeCoast.beams[1].elements[e].x1 for e in 1:nElem]

# Initialize outputs
u1 = fill([NaN], length(θRange), length(URange))
u2 = fill([NaN], length(θRange), length(URange))
u3 = fill([NaN], length(θRange), length(URange))
p1 = fill([NaN], length(θRange), length(URange))
p2 = fill([NaN], length(θRange), length(URange))
p3 = fill([NaN], length(θRange), length(URange))
M2 = fill([NaN], length(θRange), length(URange))
α = fill([NaN], length(θRange), length(URange))
cn = fill([NaN], length(θRange), length(URange))
pHinge = fill([NaN], length(θRange), length(URange))
ϕHinge = fill(NaN, length(θRange), length(URange))
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))

# Loop root pitch angle
for (j,θ) in enumerate(θRange)
    # Loop airspeed
    for (k,U) in enumerate(URange)
        println("Solving for θ = $(round(θ*180/π,digits=1)) deg, U = $U m/s")
        # Update model
        model = create_HealyLCOFFWT(wingtipTab=wingtipTab,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)
        # Set initial guess solution as the one from previous airspeed
        x0 = (k==1 || !problem[j,k-1].systemSolver.convergedFinalSolution) ? zeros(0) : problem[j,k-1].x
        # Create and solve problem
        problem[j,k] = create_SteadyProblem(model=model,systemSolver=NR,x0=x0)
        solve!(problem[j,k])
        if !(problem[j,k].systemSolver.convergedFinalSolution)
            continue
        end
        # Get outputs
        u1[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].u_n1[1],problem[j,k].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
        u2[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].u_n1[2],problem[j,k].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
        u3[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].u_n1[3],problem[j,k].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
        p1[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].p_n1[1],problem[j,k].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
        p2[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].p_n1[2],problem[j,k].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
        p3[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].p_n1[3],problem[j,k].nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)
        M2[j,k] = vcat([vcat(problem[j,k].nodalStatesOverσ[end][e].M_n1[2],problem[j,k].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
        α[j,k] = [problem[j,k].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:nElem]
        cn[j,k] = [problem[j,k].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:nElem]
        pHinge[j,k] = problem[j,k].model.hingeAxisConstraints[1].pH
        ϕHinge[j,k] = problem[j,k].model.hingeAxisConstraints[1].ϕ*180/π
    end
end

# Load reference data
aoa_25_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_2.5_exp.txt")
aoa_50_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_5.0_exp.txt")
aoa_75_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_7.5_exp.txt")
aoa_10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_10.0_exp.txt")
aoa_25_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_2.5_num.txt")
aoa_50_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_5.0_num.txt")
aoa_75_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_7.5_num.txt")
aoa_10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_10.0_num.txt")

println("Finished HealyLCOFFWTsteadyURangeAoARangeCoast.jl")