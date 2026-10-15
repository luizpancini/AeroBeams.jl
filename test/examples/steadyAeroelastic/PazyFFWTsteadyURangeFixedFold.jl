using AeroBeams

# Hinge node, fold angle [rad] and flare angle [rad]
hingeNode = 13
foldAngle = deg2rad(-90)
flareAngle = deg2rad(20)

# Root pitch angle
θ = deg2rad(7)

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:80)

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Initialize model
PazyFFWTsteadyURangeFixedFold,nElem,L,_ = create_Pazy(hinged=true,hingeNode=hingeNode,foldAngle=foldAngle,flareAngle=flareAngle,pitchAngle=θ)

# Undeformed positions and coordinates
beam = PazyFFWTsteadyURangeFixedFold.beams[1]
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
u1 = Array{Vector{Float64}}(undef,length(URange))
u2 = Array{Vector{Float64}}(undef,length(URange))
u3 = Array{Vector{Float64}}(undef,length(URange))
p1 = Array{Vector{Float64}}(undef,length(URange))
p2 = Array{Vector{Float64}}(undef,length(URange))
p3 = Array{Vector{Float64}}(undef,length(URange))
F3 = Array{Vector{Float64}}(undef,length(URange))
M2 = Array{Vector{Float64}}(undef,length(URange))
α = Array{Vector{Float64}}(undef,length(URange))
cn = Array{Vector{Float64}}(undef,length(URange))
ϕHinge = Array{Float64}(undef,length(URange))
hingeMoment = Array{Float64}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Loop airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed
    update_v_A!(PazyFFWTsteadyURangeFixedFold, v_A=[0,U,0])
    # Create and solve problem
    problem = create_SteadyProblem(model=PazyFFWTsteadyURangeFixedFold,systemSolver=systemSolver,x0=x0)
    solve!(problem)
    # Get outputs
    state = problem.loadstepHistory.states[end]
    hacState = problem.loadstepHistory.hacStates[end][1]
    aeroVars = problem.loadstepHistory.aeroVariables[end]
    u1[i] = reduce(vcat,[state[e].nodal.u1[1],
                    state[e].nodal.u2[1]] for e in 1:nElem)
    u2[i] = reduce(vcat,[state[e].nodal.u1[2],
                    state[e].nodal.u2[2]] for e in 1:nElem)                  
    u3[i] = reduce(vcat,[state[e].nodal.u1[3],
                    state[e].nodal.u2[3]] for e in 1:nElem)
    p1[i] = reduce(vcat,[state[e].nodal.p1[1],
                    state[e].nodal.p2[1]] for e in 1:nElem)
    p2[i] = reduce(vcat,[state[e].nodal.p1[2],
                    state[e].nodal.p2[2]] for e in 1:nElem)                  
    p3[i] = reduce(vcat,[state[e].nodal.p1[3],
                    state[e].nodal.p2[3]] for e in 1:nElem)
    F3[i] = reduce(vcat,[state[e].nodal.F1[3],
                    state[e].nodal.F2[3]] for e in 1:nElem)
    M2[i] = reduce(vcat,[state[e].nodal.M1[2],
                    state[e].nodal.M2[2]] for e in 1:nElem)
    α[i] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
    cn[i] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
    ϕHinge[i] = rad2deg(hacState.ϕ)
    hingeMoment[i] = hacState.hingeMoment
    # Set initial guess for next iteration
    global x0 = problem.systemArrays.x
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyURangeFixedFold"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan seen from x2-plane
plt1 = plot(xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/L, (r_n3.+u3[i])/L, aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_deformed_x2.pdf"))

# Normalized deformed wingspan seen from x1-plane
plt11 = plot(xlabel="Normalized in-plane position", ylabel="Normalized out-of-plane position", xlims=[-0.1, 0.1])
for (i,U) in enumerate(URange)
    plot!((r_n2.+u2[i])/L, (r_n3.+u3[i])/L, aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt11)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_deformed_x1.pdf"))

# Normalized deformed wingspan seen from x3-plane
plt12 = plot(xlabel="Normalized spanwise position", ylabel="Normalized in-plane position", xlims=[0, 1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/L, (r_n2.+u2[i])/L, aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt12)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_deformed_x3.pdf"))

# AoA
plt2 = plot(xlabel="Normalized spanwise position", ylabel="\$\\alpha\$ [deg]")
for (i,U) in enumerate(URange)
    plot!(x1_e/L, α[i], lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_aoa.pdf"))

# cn
plt3 = plot(xlabel="Normalized spanwise position", ylabel="\$c_n\$")
for (i,U) in enumerate(URange)
    plot!(x1_e/L, cn[i], lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_cn.pdf"))

# Hinge moment
plt5 = plot(xlabel="Airspeed [m/s]", ylabel="Hinge moment [Nm]")
plot!(URange, hingeMoment, lw=lw, c=:black, label=false)
display(plt5)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_hinge_moment.pdf"))

println("Finished PazyFFWTsteadyURangeFixedFold.jl")