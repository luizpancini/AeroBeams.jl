using AeroBeams

# Hinge node
hingeNode = 13

# Flare angle [rad]
flareAngle = deg2rad(20)

# Root pitch angle range [rad]
θRange = deg2rad.([0,3,5,7])

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:80)

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Initialize model
PazyFFWTsteadyURangeAoARangeCoast,nElem,L,_ = create_Pazy(hinged=true,hingeNode=hingeNode,flareAngle=flareAngle,pitchAngle=0)

# Undeformed positions and coordinates
beam = PazyFFWTsteadyURangeAoARangeCoast.beams[1]
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)
x1_n = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                    beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)
x1_e = [beam.elements[e].geometry.x1 for e in 1:nElem]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
F3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
M2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
α = Array{Vector{Float64}}(undef,length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(θRange),length(URange))

# Loop root pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Loop airspeed
    for (j,U) in enumerate(URange)
        println("Solving for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Update model
        model,_ = create_Pazy(hinged=true,hingeNode=hingeNode,flareAngle=flareAngle,airspeed=U,pitchAngle=θ)
        # Create and solve problem
        problem = create_SteadyProblem(model=model,systemSolver=systemSolver,x0=x0)
        solve!(problem)
        # Get outputs
        state = problem.loadstepHistory.states[end]
        hacState = problem.loadstepHistory.hacStates[end][1]
        aeroVars = problem.loadstepHistory.aeroVariables[end]
        u1[i,j] = reduce(vcat,[state[e].nodal.u1[1],
                        state[e].nodal.u2[1]] for e in 1:nElem)
        u2[i,j] = reduce(vcat,[state[e].nodal.u1[2],
                        state[e].nodal.u2[2]] for e in 1:nElem)                  
        u3[i,j] = reduce(vcat,[state[e].nodal.u1[3],
                        state[e].nodal.u2[3]] for e in 1:nElem)
        p1[i,j] = reduce(vcat,[state[e].nodal.p1[1],
                        state[e].nodal.p2[1]] for e in 1:nElem)
        p2[i,j] = reduce(vcat,[state[e].nodal.p1[2],
                        state[e].nodal.p2[2]] for e in 1:nElem)                  
        p3[i,j] = reduce(vcat,[state[e].nodal.p1[3],
                        state[e].nodal.p2[3]] for e in 1:nElem)
        F3[i,j] = reduce(vcat,[state[e].nodal.F1[3],
                        state[e].nodal.F2[3]] for e in 1:nElem)
        M2[i,j] = reduce(vcat,[state[e].nodal.M1[2],
                        state[e].nodal.M2[2]] for e in 1:nElem)
        α[i,j] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
        cn[i,j] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
        ϕHinge[i,j] = -rad2deg(hacState.ϕ)
        # Set initial guess for next iteration
        x0 = problem.systemArrays.x
    end
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyURangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
gr()

# Hinge moment
plt_HM = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]")
for (i,θ) in enumerate(θRange)
    plot!(URange, ϕHinge[i,:], lw=lw, c=colors[i], label="θ = $(round(rad2deg(θ),digits=1)) deg")
end
display(plt_HM)
savefig(string(absPath,"/PazyFFWTsteadyURangeAoARangeCoast_HM.pdf"))

println("Finished PazyFFWTsteadyURangeAoARangeCoast.jl")