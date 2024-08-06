using AeroBeams, LinearAlgebra, Plots

# Beam 
L = 2
EIy = 1e6
nElem = 20
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(EIy=EIy)],hingedNodes=[midElem+1],hingedNodesDoF=[[false,true,false]])

# BCs
q₀ = 1e0
q = (x1,t) -> ifelse.(x1.>=L/2,-q₀,0)
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
hingedBeam = create_Model(name="hingedBeam",beams=[beam],BCs=[pin,clamp])

# Create and solve the problem
problem = create_SteadyProblem(model=hingedBeam,getLinearSolution=true)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
p2 = vcat([vcat(problem.nodalStatesOverσ[end][e].p_n1[2],problem.nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Display midpoint displacement
println("Midpoint u3 = $(u3[midElem+1])")

# Plots
# ------------------------------------------------------------------------------
# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e6,save=true,savePath="/test/outputs/figures/hingedBeam/hingedBeam_deformation.pdf")
display(deformationPlot)
# u3
gr()
plt1 = plot(x1/L, u3/L, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/hingedBeam/hingedBeam_u3.pdf"))
# p2
plt2 = plot(x1/L, p2, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$p_2\$")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/hingedBeam/hingedBeam_p2.pdf"))
# F3
plt3 = plot(x1/L, F3/1e3, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [kN]")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/hingedBeam/hingedBeam_F3.pdf"))
# M2
plt4 = plot(x1/L, M2/1e3, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [kN.m]")
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/hingedBeam/hingedBeam_M2.pdf"))

println("Finished hingedBeam.jl")