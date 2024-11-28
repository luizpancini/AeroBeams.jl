using AeroBeams

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 12
foldAngle = -π/2*1/2
flareAngle = 20*π/180

# Spring stiffness
kSpring = 0e6

# Root pitch angle
θ = 5*pi/180

# Gravity
g = 9.80665

# Airspeed range
URange = collect(1:1:40)

# Initialize model
pazyFFWT = create_PazyFFWT(hingeNode=hingeNode,foldAngle=foldAngle,flareAngle=flareAngle,kSpring=kSpring,p0=[0;0;θ],g=g)

# System solver
σ0 = 1
maxIter = 50
relTol = 1e-5
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Fixed properties of the model
elemNodes = vcat([vcat(pazyFFWT.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [pazyFFWT.r_n[n][1] for n in elemNodes]
r_n2 = [pazyFFWT.r_n[n][2] for n in elemNodes]
r_n3 = [pazyFFWT.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(pazyFFWT.beams[1].elements[e].x1_n1,pazyFFWT.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [pazyFFWT.beams[1].elements[e].x1 for e in 1:15]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(URange))
u3 = Array{Vector{Float64}}(undef,length(URange))
p1 = Array{Vector{Float64}}(undef,length(URange))
p2 = Array{Vector{Float64}}(undef,length(URange))
M2 = Array{Vector{Float64}}(undef,length(URange))
α = Array{Vector{Float64}}(undef,length(URange))
cn = Array{Vector{Float64}}(undef,length(URange))
problem = Array{SteadyProblem}(undef,length(URange))

# Loop airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed
    set_motion_basis_A!(model=pazyFFWT,v_A=[0;U;0])
    # Set initial guess solution as the one from previous airspeed
    x0 = (i==1) ? zeros(0) : problem[i-1].x
    # Create and solve problem
    problem[i] = create_SteadyProblem(model=pazyFFWT,systemSolver=NR,x0=x0)
    solve!(problem[i])
    # Get outputs
    u1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...)
    u3[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
    p1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].p_n1[1],problem[i].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
    p2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].p_n1[2],problem[i].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
    M2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].M_n1[2],problem[i].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)
    α[i] = [problem[i].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
    cn[i] = [problem[i].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]
end

println("Finished PazyFFWTsteadyURange.jl")