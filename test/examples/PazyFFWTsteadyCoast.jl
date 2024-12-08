using AeroBeams

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 13
flareAngle = 10*π/180

# Spring stiffness
kSpring = 0e-4

# Root pitch angle
θ = 7*pi/180

# Gravity
g = 9.8

# Airspeed
U = 20

# Pazy wing with flared folding tip
PazyFFWTsteadyCoast = create_PazyFFWT(hingeNode=hingeNode,flareAngle=flareAngle,kSpring=kSpring,airspeed=U,p0=[0;0;θ],g=g)

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 100
relTol = 1e-9
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,minimumLoadFactorStep=σstep,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

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
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
p1 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[1],problem.nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
p3 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[3],problem.nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)
α = [problem.aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
cn = [problem.aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]

# Get rotation at the hinge node
ΔϕHinge = problem.model.hingeAxisConstraints[1].Δϕ*180/π

# Display results
println("Coast angle = $ΔϕHinge deg")

println("Finished PazyFFWTsteadyCoast.jl")