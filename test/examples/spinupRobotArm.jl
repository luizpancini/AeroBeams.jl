using AeroBeams, LinearAlgebra

# Select to promote rotation through basis A (does not work as expected otherwise)
rotateBasisA = true

# Rotation variables
τ = 15.0
θ = t -> ifelse.(t.<=τ, 6/15*(t.^2/2 .+ (15/(2*π))^2 * (cos.(2*π*t/τ).-1)), 6*t.-45)
θdot = t -> ifelse.(t.<=τ, 2*t/5 .- 3/π*sin.(2*pi*t/τ), 6)
p = t -> 4*tan.(θ(t)/4)
pdot = t -> sec.(θ(t)/4).^2 .* θdot(t)

# Beam
L = 10
EA,GA,GJ,EI = 2.8e7,1e7,1.4e4,1.4e4
ρA,ρI = 1.2,6e-4
nElements = 3
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElements,C=[stiffnessMatrix],I=[inertiaMatrix])
if !rotateBasisA
    add_initial_displacements_and_velocities_to_beam!(beam,conditionTypes=["pdot0_of_x1"],conditionFuns=[(x1)->[0; 0; pdot(0)]])
end

# BCs
driver = create_BC(name="driver",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,t->p(t)])
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A"],values=[0,0,0,0,0])

# Model
if rotateBasisA
    spinupRobotArm = create_Model(name="spinupRobotArm",beams=[beam],BCs=[pin],ω_A=t->[0;0;-θdot(t)])
else
    spinupRobotArm = create_Model(name="spinupRobotArm",beams=[beam],BCs=[driver])
end

# Time variables
tf = 20
Δt = 4e-2
time = unique(vcat(collect(0:Δt/2:tf/2),collect(tf/2:Δt:tf)))

# Create and solve the problem
problem = create_DynamicProblem(model=spinupRobotArm,timeVector=time)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u1_tip = [problem.nodalStatesOverTime[i][nElements].u_n2[1] for i in 1:length(t)]
u2_tip = [problem.nodalStatesOverTime[i][nElements].u_n2[2] for i in 1:length(t)]
θ3_root = [problem.nodalStatesOverTime[i][1].θ_n1 for i in 1:length(t)]

# Rigid beam's tip displacements
u1_tip_rigid = L*(cos.(θ(t)).-1)
u2_tip_rigid = L*sin.(θ(t))

println("Finished spinupRobotArm.jl")