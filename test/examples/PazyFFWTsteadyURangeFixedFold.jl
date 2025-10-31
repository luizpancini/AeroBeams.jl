using AeroBeams, LinearAlgebra

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 13
foldAngle = -30*π/180
flareAngle = 20*π/180

# Spring stiffness
kSpring = 1e-4
kIPBendingHinge = 1e-1

# Root pitch angle
θ = 7*pi/180

# Airspeed range
URange = collect(1:1:60)

# Solution method for constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = false

# Flag to use previous solution as initial guess
usePreviousSolution = true

# Initialize model
PazyFFWTsteadyURangeFixedFold = create_PazyFFWT(hingeNode=hingeNode,foldAngle=foldAngle,flareAngle=flareAngle,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,pitchAngle=θ)

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-6
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Fixed properties of the model
elemNodes = vcat([vcat(PazyFFWTsteadyURangeFixedFold.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [PazyFFWTsteadyURangeFixedFold.r_n[n][1] for n in elemNodes]
r_n2 = [PazyFFWTsteadyURangeFixedFold.r_n[n][2] for n in elemNodes]
r_n3 = [PazyFFWTsteadyURangeFixedFold.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(PazyFFWTsteadyURangeFixedFold.beams[1].elements[e].x1_n1,PazyFFWTsteadyURangeFixedFold.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [PazyFFWTsteadyURangeFixedFold.beams[1].elements[e].x1 for e in 1:15]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(URange))
u2 = Array{Vector{Float64}}(undef,length(URange))
u3 = Array{Vector{Float64}}(undef,length(URange))
p1 = Array{Vector{Float64}}(undef,length(URange))
p2 = Array{Vector{Float64}}(undef,length(URange))
p3 = Array{Vector{Float64}}(undef,length(URange))
M2 = Array{Vector{Float64}}(undef,length(URange))
α = Array{Vector{Float64}}(undef,length(URange))
cn = Array{Vector{Float64}}(undef,length(URange))
pHinge = Array{Vector{Float64}}(undef,length(URange))
ϕHinge = Array{Float64}(undef,length(URange))
hingeMoment = Array{Float64}(undef,length(URange))
problem = Array{SteadyProblem}(undef,length(URange))

# Loop airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed
    model = create_PazyFFWT(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,hingeNode=hingeNode,foldAngle=foldAngle,flareAngle=flareAngle,kSpring=kSpring,kIPBendingHinge=kIPBendingHinge,pitchAngle=θ,airspeed=U)
    # Set initial guess solution as the one from previous airspeed
    if usePreviousSolution
        x0 = (i==1) ? zeros(0) : problem[i-1].x
    else
        x0 = zeros(0)
    end
    # Create and solve problem
    problem[i] = create_SteadyProblem(model=model,systemSolver=NR,x0=x0)
    solve!(problem[i])
    # Get outputs
    u1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...)
    u2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[2],problem[i].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:15]...)
    u3[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...)
    p1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].p_n1[1],problem[i].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...)
    p2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].p_n1[2],problem[i].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...)
    p3[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].p_n1[3],problem[i].nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...)
    M2[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].M_n1[2],problem[i].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...)
    α[i] = [problem[i].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15]
    cn[i] = [problem[i].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15]
    pHinge[i] = problem[i].model.hingeAxisConstraints[1].pH
    ϕHinge[i] = problem[i].model.hingeAxisConstraints[1].ϕ*180/π
    hingeMoment[i] = problem[i].model.hingeAxisConstraints[1].hingeMoment
end

println("Finished PazyFFWTsteadyURangeFixedFold.jl")