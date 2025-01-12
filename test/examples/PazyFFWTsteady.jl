using AeroBeams

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 13
foldAngle = -π/2*2/2
flareAngle = π*20/180

# Spring stiffness
kSpring = 0

# Root pitch angle
θ = 7*pi/180

# Airspeed
U = 50

# Gravity
g = 9.8

# Pazy wing with flared folding tip
pazyFFWT = create_PazyFFWT(hingeNode=hingeNode,foldAngle=foldAngle,flareAngle=flareAngle,kSpring=kSpring,airspeed=U,p0=[0;0;θ],g=g)

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 100
relTol = 1e-6
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter,relativeTolerance=relTol)

# Create and solve problem
problem = create_SteadyProblem(model=pazyFFWT,systemSolver=NR)
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

# Compute rotations, displacement and balance moment at the hinge node
Δp_WMHingeNode = [p1[2*hingeNode]-p1[2*hingeNode-2]; p2[2*hingeNode]-p2[2*hingeNode-2]; p3[2*hingeNode]-p3[2*hingeNode-2]]
R,_ = rotation_tensor_WM(Δp_WMHingeNode)
yaw,pitch,roll = ypr_from_rotation_tensor(R,assumeNullYawInSingularity=true)
hingeAngles = [yaw,pitch,roll]*180/π
u3HingeNode = u3[2*hingeNode-1]
hingeBalanceM2 = -problem.model.rotationConstraints[1].balanceMoment
hingeBalanceM1 = -problem.model.rotationConstraints[2].balanceMoment
hingeBalanceM3 = -problem.model.rotationConstraints[3].balanceMoment
hingeBalanceM = sqrt(hingeBalanceM1^2+hingeBalanceM2^2+hingeBalanceM3^2)

# Display results
println("Twist rotation at hinge node = $(hingeAngles[3]) deg")
println("OOP rotation at hinge node = $(hingeAngles[2]) deg")
println("IP rotation at hinge node = $(hingeAngles[1]) deg")
println("OOP displacement at hinge node = $u3HingeNode")
println("Moment necessary to impose the hinge angle = $hingeBalanceM")

println("Finished PazyFFWTsteady.jl")