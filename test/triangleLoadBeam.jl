using AeroBeams, LinearAlgebra, Plots

# Beam
L = 1.0
EIy = 1e7
∞ = 1e12
stiffnessMatrix = diagm([∞,∞,∞,∞,EIy,∞])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
q₀ = 1e6
q = (x1,t) -> q₀*x1
add_loads_to_beam!(beam,loadTypes=["f_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q(x1,t)]])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
roller = create_BC(name="roller",beam=beam,node=1,types=["u3A"],values=[0])

# Model
triangleLoadBeam = create_Model(name="triangleLoadBeam",beams=[beam],BCs=[clamp,roller])

# Create and solve the problem
problem = create_SteadyProblem(model=triangleLoadBeam,getLinearSolution=true)
solve!(problem)

# Get nodal arclength positions, displacements and forces over the beams
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Plots
# ------------------------------------------------------------------------------
# Deformed shape
deformationPlot = plot_steady_deformation(problem,scale=1e3,save=true,savePath="/test/outputs/figures/triangleLoadBeam/triangleLoadBeam_deformation.pdf")
display(deformationPlot)
# Displacement
gr()
plt1 = plot(x1/L, u3, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3\$ [m]")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/triangleLoadBeam/triangleLoadBeam_u3.pdf"))
# Force
plt2 = plot(x1/L, F3, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/triangleLoadBeam/triangleLoadBeam_F3.pdf"))
# Moment
plt3 = plot(x1/L, M2, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/triangleLoadBeam/triangleLoadBeam_M2.pdf"))

println("Finished triangleLoadBeam.jl")