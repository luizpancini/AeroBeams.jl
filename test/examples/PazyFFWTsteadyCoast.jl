using AeroBeams

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 12
flareAngle = 20*π/180

# Spring stiffness
kSpring = 1e-4
kIPBendingHinge = 1e0

# Root pitch angle
θ = 7*π/180

# Airspeed
U = 35

# Solution method for constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = false

# Pazy wing with flared folding tip
PazyFFWTsteadyCoast = create_PazyFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,hingeNode=hingeNode,flareAngle=flareAngle,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ)

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=PazyFFWTsteadyCoast,systemSolver=NR)
solve!(problem)

# Get outputs
elemNodes = vcat([vcat(problem.model.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [problem.model.r_n[n][1] for n in elemNodes]
r_n2 = [problem.model.r_n[n][2] for n in elemNodes]
r_n3 = [problem.model.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [problem.model.beams[1].elements[e].x1 for e in 1:15]
u1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...)
u2 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[2],problem.nodalStatesOverσ[end][e].u_n2[2]) for e in 1:15]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
p3 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[3],problem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)
α = [problem.aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
cn = [problem.aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]

# Get rotation at the hinge node
ϕHinge = problem.model.hingeAxisConstraints[1].ϕ*180/π

# Display results
println("Coast angle = $ϕHinge deg")

println("Finished PazyFFWTsteadyCoast.jl")