using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam 
R,θ = 2.54,120*π/180
L = R*θ
A,Iy = 4.05e-4,13.1e-8
E = 70.4e9
∞ = 1e12
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 80
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-θ/2;0],k=[0;1/R;0])

# BCs
λ = 8.9
q = -λ*EIy/R^2
add_loads_to_beam!(beam,loadTypes=["ff_A_of_x1t"],loadFuns=[(x1,t)->[0; 0; q]])
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
archUnderFollowerPressure = Model(name="archUnderFollowerPressure",beams=[beam],BCs=[clamp1,clamp2])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = SteadyProblem(model=archUnderFollowerPressure,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
mid_u3 = [problem.nodalStatesOverσ[i][div(nElem,2)].u_n2[3] for i in 1:length(σVector)]

# Plot normalized displacements over load steps
plt1 = plot(-mid_u3/R, σVector*λ, color=:black, lw=2, xlabel="Midpoint \$-u_3/R\$", ylabel="\$\\lambda\$", label=false)
display(plt1)

println("Finished archUnderFollowerPressure.jl")