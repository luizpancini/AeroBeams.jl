using AeroBeams

# Hinge node and flare angle [rad]
hingeNode = 13
flareAngle = 20*π/180

# Spring stiffness
kSpring = 1e-4

# Gravity
g = 9.8

# Root pitch angle range
θRange = π/180*[0,3,5,7]

# Airspeed range
URange = collect(1:1:80)

# Initialize model
PazyFFWTsteadyURangeAoARangeCoast = create_PazyFFWT(hingeNode=hingeNode,flareAngle=flareAngle,kSpring=kSpring,pitchAngle=θ,g=g)

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-5
ΔλRelaxFactor = 1
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol,ΔλRelaxFactor=ΔλRelaxFactor)

# Fixed properties of the model
elemNodes = vcat([vcat(PazyFFWTsteadyURangeAoARangeCoast.elements[e].nodesGlobalID) for e in 1:15]...)
r_n1 = [PazyFFWTsteadyURangeAoARangeCoast.r_n[n][1] for n in elemNodes]
r_n2 = [PazyFFWTsteadyURangeAoARangeCoast.r_n[n][2] for n in elemNodes]
r_n3 = [PazyFFWTsteadyURangeAoARangeCoast.r_n[n][3] for n in elemNodes]
x1_n = vcat([vcat(PazyFFWTsteadyURangeAoARangeCoast.beams[1].elements[e].x1_n1,PazyFFWTsteadyURangeAoARangeCoast.beams[1].elements[e].x1_n2) for e in 1:15]...)
x1_e = [PazyFFWTsteadyURangeAoARangeCoast.beams[1].elements[e].x1 for e in 1:15]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
M2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
α = Array{Vector{Float64}}(undef,length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(URange))
pHinge = Array{Vector{Float64}}(undef,length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(θRange),length(URange))
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))

# Loop root pitch angle
for (i,θ) in enumerate(θRange)
    # Loop airspeed
    for (j,U) in enumerate(URange)
        println("Solving for θ = $(θ*180/π) deg, U = $U m/s")
        # Update model
        model = create_PazyFFWT(hingeNode=hingeNode,flareAngle=flareAngle,kSpring=kSpring,airspeed=U,pitchAngle=θ,g=g)
        # Set initial guess solution as the one from previous airspeed
        x0 = (j==1) ? zeros(0) : problem[i,j-1].x
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,systemSolver=NR)
        solve!(problem[i,j])
        # Get TF for converged solution
        converged = problem[i,j].systemSolver.convergedFinalSolution
        # Get outputs
        u1[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[1],problem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:15]...) : fill(NaN64,30)
        u2[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[2],problem[i,j].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:15]...) : fill(NaN64,30)
        u3[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[3],problem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:15]...) : fill(NaN64,30)
        p1[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].p_n1[1],problem[i,j].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:15]...) : fill(NaN64,30)
        p2[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].p_n1[2],problem[i,j].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:15]...) : fill(NaN64,30)
        p3[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].p_n1[3],problem[i,j].nodalStatesOverσ[end][e].p_n2[3]) for e in 1:15]...) : fill(NaN64,30)
        M2[i,j] = converged ? vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].M_n1[2],problem[i,j].nodalStatesOverσ[end][e].M_n2[2]) for e in 1:15]...) : fill(NaN64,30)
        α[i,j] = converged ? [problem[i,j].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:15] : fill(NaN64,15)
        cn[i,j] = converged ? [problem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:15] : fill(NaN64,15)
        pHinge[i,j] = converged ? problem[i,j].model.hingeAxisConstraints[1].pH : fill(NaN64,3)
        ϕHinge[i,j] = converged ? problem[i,j].model.hingeAxisConstraints[1].ϕ*180/π : NaN64
    end
end

# Reverse sign for ϕHinge at U=0, if applicable
for (i,θ) in enumerate(θRange)
    ϕHinge[i,1] = ϕHinge[i,1] < 0 ? -ϕHinge[i,1] : ϕHinge[i,1]
end    

println("Finished PazyFFWTsteadyURangeAoARangeCoast.jl")