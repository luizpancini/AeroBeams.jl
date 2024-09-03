# # Arch under follower pressure
# This example simulates the static response of an arch subjected to a normal (follower) pressure. Let's begin by loading the package.

using AeroBeams

# ### Beam 
# The first step is to create a beam. The arch has radius `R` and spans over and angle `θ`. The curvature of the beam is thus `1/R`, and the total length is `L=Rθ`. We define the beam orientation such that the arch spans from an angle `-θ/2` to `θ/2` about a vertical line. This is done by specifying the rotation parameters from basis `A` to basis `b`, `p0`, with the Euler parameters sequence 3-2-1: the angle of rotation about the second axis is `-θ/2`. The cross-section has area `A` and bending moment of inertia `Iy`. The elastic modulus of the material is `E`. We discretize the beam into `nElem` finite elements.
R,θ = 2.54,120*π/180
k2 = 1/R
L = R*θ
A,Iy = 4.05e-4,13.1e-8
E = 70.4e9
EA,EIy = E*A,E*Iy
nElem = 80
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=1e12,EA=EA,EIy=EIy)],rotationParametrization="E321",p0=[0;-θ/2;0],k=[0;k2;0]);

# BCs
λ = 11
q = -λ*EIy/R^2
add_loads_to_beam!(beam,loadTypes=["ff_b_of_x1t"],loadFuns=[(x1,t)->[0; 0; q]])
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
archUnderFollowerPressure = create_Model(name="archUnderFollowerPressure",beams=[beam],BCs=[clamp1,clamp2])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=archUnderFollowerPressure,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
mid_u3 = [problem.nodalStatesOverσ[i][div(nElem,2)].u_n2[3] for i in 1:length(σVector)]

println("Finished archUnderFollowerPressure.jl")