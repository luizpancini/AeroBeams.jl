using AeroBeams, LinearAlgebra

# Beam
L,b,H = 479e-3,50.8e-3,0.45e-3
A,Iy,Iz = b*H,b*H^3/12,H*b^3/12
J = Iy+Iz
Ksy = Ksz = 5/6
E,ν,ρ = 127e9,0.36,4.43e3
G = E/(2*(1+ν))
nElem = 60
stiffnessMatrix = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J,E*Iy,E*Iz])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*(Iy+Iz),ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,-π/2,0])

# BCs
ω = 32*(2*π)
V = 0.3414
A = V/ω
T = 2*π/ω
u₃b = t -> A*sin.(ω*t)
shaker = create_BC(name="shaker",beam=beam,node=1,types=["u1b","u2b","u3b","p1b","p2b","p3b"],values=[0,0,t->u₃b(t),0,0,0])

# Model
rootExcitationBeam2 = create_Model(name="rootExcitationBeam2",beams=[beam],BCs=[shaker],gravityVector=[0,0,-9.80665])

# Time variables
cycles = 5
tf = cycles*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=rootExcitationBeam2,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
x1 = [v.x1 for v in rootExcitationBeam2.elements]
u3b_root = [problem.nodalStatesOverTime[i][1].u_n1_b[3] for i in 1:length(t)]
u3b_tip = [problem.nodalStatesOverTime[i][nElem].u_n2_b[3] for i in 1:length(t)]
V3_root = [problem.elementalStatesOverTime[i][1].V[3] for i in 1:length(t)]
V3_tip = [problem.elementalStatesOverTime[i][nElem].V[3] for i in 1:length(t)]
V3 = Vector{Vector{Float64}}()
for i in 1:length(t)
    push!(V3,[problem.elementalStatesOverTime[i][e].V[3] for e = 1:nElem])
end

println("Finished rootExcitationBeam2.jl")