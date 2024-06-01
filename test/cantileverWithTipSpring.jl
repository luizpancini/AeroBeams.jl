using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam
L = 1
EIy = 1e4
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(EIy=EIy)])

# Spring
ku = [0; 0; 1e4]
spring = create_Spring(elementID=nElem,localNode=2,ku=ku)
add_springs_to_beam!(beam,springs=[spring])

# BCs
F = 1e3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["F3A"],values=[F])

# Model
cantileverWithTipSpring = create_Model(name="cantileverWithTipSpring",beams=[beam],BCs=[clamp,tipForce])

# Create and solve the problem
problem = create_SteadyProblem(model=cantileverWithTipSpring)
solve!(problem)

# Get nodal arclength positions, displacement and loads over the beam
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)

# Tip displacement
println("Tip disp = $(u3[end])")

# Plots
plt1 = plot(x1/L, u3/L, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3/L\$")
display(plt1)
plt2 = plot(x1/L, F3, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
display(plt2)
plt3 = plot(x1/L, M2, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
display(plt3)

println("Finished cantileverWithTipSpring.jl")