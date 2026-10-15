using AeroBeams

# Beam 
L = 1
EIy = 1
nElem = 40
midElem = div(nElem,2)
hingedNode = midElem+1
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[hingedNode],hingedNodesDoF=[[false,true,false]])

# Spring
k = 1e0
springOOP = create_Spring(kp=k)
add_spring_to_beams!(beams=[beam,beam],spring=springOOP,nodesLocalIDs=[hingedNode-1,hingedNode+1],p=[π/2,0,0])

# BCs (remove inboard moment to get deformed inner shape)
M = π/2
inboardMoment = create_BC(name="moment",beam=beam,node=hingedNode-1,conditions=Dict(:M2A=>-M))
outboardMoment = create_BC(name="moment",beam=beam,node=hingedNode+1,conditions=Dict(:M2A=>M))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
drivenSpringedHingedBeam2 = create_Model(name="drivenSpringedHingedBeam2",beams=[beam],BCs=[clamp,inboardMoment,outboardMoment])

# System solver
systemSolver = create_NLsolver(displayStatus=true,method=:newton,linesearch=true)

# Create and solve problem
problem = create_SteadyProblem(model=drivenSpringedHingedBeam2,systemSolver=systemSolver)
solve!(problem)

# Undeformed positions and coordinates
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)
x1_n = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                    beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)

# Solution outputs
state = problem.loadstepHistory.states[end]
u1 = reduce(vcat,[state[e].nodal.u1[1],
                  state[e].nodal.u2[1]] for e in 1:nElem)
u2 = reduce(vcat,[state[e].nodal.u1[2],
                  state[e].nodal.u2[2]] for e in 1:nElem)                  
u3 = reduce(vcat,[state[e].nodal.u1[3],
                  state[e].nodal.u2[3]] for e in 1:nElem)
p1 = reduce(vcat,[state[e].nodal.p1[1],
                  state[e].nodal.p2[1]] for e in 1:nElem)
p2 = reduce(vcat,[state[e].nodal.p1[2],
                  state[e].nodal.p2[2]] for e in 1:nElem)                  
p3 = reduce(vcat,[state[e].nodal.p1[3],
                  state[e].nodal.p2[3]] for e in 1:nElem)
F3 = reduce(vcat,[state[e].nodal.F1[3],
                  state[e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[state[e].nodal.M1[2],
                  state[e].nodal.M2[2]] for e in 1:nElem)

## Plots
# Set paths
relPath = "/test/outputs/figures/drivenSpringedHingedBeam2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Plot rotations and bending moments
plt_p2 = plot(x1_n,p2,ylabel="\$p_2\$",lw=lw,label=false)
plt_M2 = plot(x1_n,M2,xlabel="\$x_1\$",ylabel="\$M_2\$ [N.m]",lw=lw,label=false)
plt_pM = plot(plt_p2, plt_M2, layout=(2,1))
display(plt_pM)
savefig(string(absPath,"/drivenSpringedHingedBeam2_pM.pdf"))

# Plot 2D position (x1-x3 plane)
plt_pos_x1x3 = plot(r_n1+u1, r_n3+u3, xlabel="\$x_1\$", ylabel="\$x_3\$", lw=lw, label=false, aspect_ratio=:equal)
display(plt_pos_x1x3)
savefig(string(absPath,"/drivenSpringedHingedBeam2_pos_x1x3.pdf"))

println("Finished drivenSpringedHingedBeam2.jl")