using AeroBeams, LinearAlgebra, DelimitedFiles
 
# Beam
L = 0.5
EA,GA,GJ,EIy,EIz = 2e10*1e-6,1e12,1e12,1e12,1e12
ρA,ρI = 8e9*1e-6,0
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# Time variables
tf = 1.0
Δt = 2e-3

# BCs
F₀ = 1.0
F = t -> ifelse.(t.>=Δt, F₀, 0.0)
axialForce = create_BC(name="axialForce",beam=beam,node=nElem+1,types=["F1A"],values=[t->F(t)])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
axialTractionCantilever = create_Model(name="axialTractionCantilever",beams=[beam],BCs=[clamp,axialForce])

# Create and solve the problem
problem = create_DynamicProblem(model=axialTractionCantilever,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
x1 = [axialTractionCantilever.r_n[i][1] for i in 1:nElem+1]
ind_t_10 = length(problem.timeVector)
ind_t_08 = round(Int,0.8*ind_t_10)
u1_08 = [problem.nodalStatesOverTime[ind_t_08][i].u_n2[1] for i in 1:nElem]
u1_10 = [problem.nodalStatesOverTime[ind_t_10][i].u_n2[1] for i in 1:nElem]
insert!(u1_08, 1, 0)
insert!(u1_10, 1, 0)

# Load reference solution (by Reddy)
x1_ref = collect(LinRange(0,L,21))
u1_08_ref = vec(readdlm(joinpath(dirname(@__DIR__), "referenceData", "axialTractionCantilever", "t08.txt")))
u1_10_ref = vec(readdlm(joinpath(dirname(@__DIR__), "referenceData", "axialTractionCantilever", "t08.txt")))

println("Finished axialTractionCantilever.jl")