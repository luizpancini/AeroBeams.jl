using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Geometric properties
L = 37.5
b,H = 1,0.117
A,Iy,Iz = b*H,b*H^3/12,H*b^3/12
J = Is = Iy + Iz
Ksy,Ksz,Kt = 5/6,1/14.625,1/65.852
r0 = [2.5; 0.0; 0.0]

# Material properties
E11,E22,E33 = 2.059e7,1.42e6,1.42e6
G12,G13,G23 = 8.9e5,8.9e5,8e5
ν12,ν13,ν23 = 0.42,0.42,0.54
ρ = 1.44e-4

# Discretization variables
nElem = 30

# Range of angular velocity [rad/s]
ωRPM = collect(0:50:750)
ωRange = 2*π/60*ωRPM

# Number of modes
nModes = 6

# Beams
complianceMatrix = zeros(6,6)
complianceMatrix[1:3,1:3] .= [1/(E11*A) 0 0;
                              0 1/(G12*A*Ksy) 0;
                              0 0 1/(G13*A*Ksz)]
complianceMatrix[4:6,4:6] .= diagm([1/(G23*J*Kt); 1/(E11*Iy); 1/(E11*Iz)])
stiffnessMatrix = inv(complianceMatrix)
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Create model
straightRotor = create_Model(name="straightRotor",beams=[beam],BCs=[clamp],initialPosition=r0,units=create_UnitsSystem(length="in",force="lbf",frequency="Hz"))

# Initialize outputs
numFreqs = Vector{Vector{Float64}}(undef,length(ωRange))

# Loop over angular velocities
for (i,ω) in enumerate(ωRange)
    # Set angular velocity of basis A and update model
    straightRotor.ω_A = [0; 0; ω]
    update_model!(straightRotor)
    # Create and solve eigenproblem
    global problem = create_EigenProblem(model=straightRotor,nModes=nModes)
    solve!(problem)
    # Get outputs
    numFreqs[i] = problem.frequenciesOscillatory
end

# Load experimental values from Epps & Chandra (1996)
ωRPMExp = [0,250,500,750]
expFreqs = [5.2374e+00   7.4239e+00   1.1018e+01   1.5820e+01;
            3.2999e+01   3.4582e+01   4.0188e+01   4.7001e+01;
            9.0439e+01   9.1112e+01   9.7733e+01   1.0558e+02;
            1.7237e+02   1.7532e+02   1.7967e+02   1.8874e+02]
            
println("Finished straightRotor.jl")