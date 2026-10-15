using AeroBeams

# Spring stiffness range
kRange = vcat(5e-3,1e-2,1e-1,5e-1,1,10)

# Beam 
L = 1
EIy = 1
nElem = 20
hingeNode = div(nElem,2)+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingeNode],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = -1
q = (x1,t) -> q₀
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t => (x1,t)->[0; 0; q(x1,t)]))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
clampedHingedBeamSpringRange = create_Model(name="clampedHingedBeamSpringRange",beams=[beam],BCs=[clamp])

# Initialize outputs
u1 = Vector{Vector{Float64}}(undef,length(kRange))
u3 = Vector{Vector{Float64}}(undef,length(kRange))
p2 = Vector{Vector{Float64}}(undef,length(kRange))
F3 = Vector{Vector{Float64}}(undef,length(kRange))
M2 = Vector{Vector{Float64}}(undef,length(kRange))
springs = Vector{Spring}(undef,length(kRange))
hingeAngle = Vector{Float64}(undef,length(kRange))
springMoment = Vector{Float64}(undef,length(kRange))
problem = Vector{SteadyProblem}(undef,length(kRange))

# Loop spring stiffness
for (i,k) in enumerate(kRange)
    # Show progress
    println("Solving for k = $k")
    # Remove springs from beam
    remove_all_springs_from_beams!(beams=[beam])
    # Add spring with current stiffness to beam
    springs[i] = create_Spring(kp=k)
    add_spring_to_beams!(clampedHingedBeamSpringRange,beams=[beam,beam],spring=springs[i],nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
    # Set initial guess solution as previous known solution
    x0 = i==1 ? 
         zeros(0) : 
         problem[i-1].systemArrays.x
    # Create and solve the problem
    problem[i] = create_SteadyProblem(model=clampedHingedBeamSpringRange,x0=x0)
    solve!(problem[i])
    # Get outputs
    histLast = problem[i].loadstepHistory.states[end]
    u1[i] = reduce(vcat,
            [histLast[e].nodal.u1[1],
             histLast[e].nodal.u2[1]] for e in 1:nElem)
    u3[i] = reduce(vcat,
            [histLast[e].nodal.u1[3],
             histLast[e].nodal.u2[3]] for e in 1:nElem)
    p2[i] = reduce(vcat,
            [histLast[e].nodal.p1[2],
             histLast[e].nodal.p2[2]] for e in 1:nElem)
    F3[i] = reduce(vcat,
            [histLast[e].nodal.F1[3],
             histLast[e].nodal.F2[3]] for e in 1:nElem)
    M2[i] = reduce(vcat,
            [histLast[e].nodal.M1[2],
             histLast[e].nodal.M2[2]] for e in 1:nElem)
    hingeAngle[i] = rad2deg(rotation_angle_limited([0; p2[i][2*hingeNode-1]-p2[i][2*hingeNode-2]; 0]))
    springMoment[i] = norm(springs[i].state.Ms)
end

# Get nodal arclength and coordinates
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)

## Plots
# Set paths
relPath = "/test/outputs/figures/clampedHingedBeamSpringRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, kRange)
lw = 2
gr()

# Deformed shape
plt_u = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
for (i,k) in enumerate(kRange)
    plot!((r_n1.+u1[i])/x1[end], (r_n3.+u3[i])/x1[end], lz=k, c=colors, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt_u)
savefig(string(absPath,"/clampedHingedBeamSpringRange_u.pdf"))

# F3
plt_F3 = plot(xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
for (i,k) in enumerate(kRange)
    plot!(x1/L, F3[i], lz=k, c=colors, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt_F3)
savefig(string(absPath,"/clampedHingedBeamSpringRange_F3.pdf"))

# M2
plt_M2 = plot(xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
for (i,k) in enumerate(kRange)
    plot!(x1/L, M2[i], lz=k, c=colors, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt_M2)
savefig(string(absPath,"/clampedHingedBeamSpringRange_M2.pdf"))

# Spring moment
maximumTheoreticalSpringMoment = abs(q₀)*(L/2)*(L/2)/2
plt_Ms = plot(xlabel="\$log(k)\$", ylabel="\$M_s\$ [N.m]", xlims=log10.([minimum(kRange), maximum(kRange)]), ylims=[0,maximumTheoreticalSpringMoment], yticks=0:0.025:maximumTheoreticalSpringMoment)
plot!(log10.(kRange),springMoment, lw=lw, marker=:circle, label=false)
display(plt_Ms)
savefig(string(absPath,"/clampedHingedBeamSpringRange_Ms.pdf"))

# Hinge angle
plt_phi = plot(xlabel="\$log(k)\$", ylabel="Hinge angle [deg]", xlims=log10.([minimum(kRange), maximum(kRange)]), ylims=[0,90], yticks=[-90,-45,0,45,90])
plot!(log10.(kRange),hingeAngle, lw=lw, marker=:circle, label=false)
display(plt_phi)
savefig(string(absPath,"/clampedHingedBeamSpringRange_phi.pdf"))

println("Finished clampedHingedBeamSpringRange.jl")