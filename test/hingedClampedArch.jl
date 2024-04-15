using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam 
R,θ= 100,215*π/180
L = R*θ
A,Iy = 2,2/3
E = 1.5e6
∞ = 1e10
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 80
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-θ/2;0],k=[0;1/R;0])

# BCs
λ = 8.9
elemForce = div(nElem,2)
hinge = create_BC(name="hinge",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
force = create_BC(name="force",beam=beam,node=elemForce+1,types=["F3A"],values=[-λ*E*Iy/R^2])

# Model
hingedClampedArch = Model(name="hingedClampedArch",beams=[beam],BCs=[hinge,clamp,force])

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = SteadyProblem(model=hingedClampedArch,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
u1_atForce = [problem.nodalStatesOverσ[i][elemForce].u_n2[1] for i in 1:length(σVector)]
u3_atForce = [problem.nodalStatesOverσ[i][elemForce].u_n2[3] for i in 1:length(σVector)]

# Plot normalized displacements over load steps
x = [-u1_atForce/R, -u3_atForce/R]
labels = ["\$-u_1/R\$" "\$-u_3/R\$"]
colors = [:blue,:orange]
plt1 = plot()
plot!(x, σVector*λ, linewidth=2, label=false, xlabel="\$-u_1/L, -u_3/L,\$", ylabel="\$\\lambda\$", title="Displacements at point of force application")
halfNσ = round(Int,length(σVector)/2)
annotate!(x[1][halfNσ], σVector[halfNσ]*λ, text(labels[1], :bottom, :right, colors[1]))
annotate!(x[2][halfNσ], σVector[halfNσ]*λ, text(labels[2], :top, :left, colors[2]))
display(plt1)

println("Finished hingedClampedArch.jl")