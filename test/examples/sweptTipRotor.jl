using AeroBeams, LinearAlgebra

# Geometric properties
L1,L2 = 31.5,6.0
b,H = 1,0.063
A,Iy,Iz = b*H,b*H^3/12,H*b^3/12
J = Is = Iy + Iz
Ksy,Ksz,Kt = 5/6,1/14.625,1/65.852
r0 = [2.5; 0.0; 0.0]

# Material properties
E = 1.06e7
G = E/(2*(1+0.325))
ρ = 2.51e-4

# Discretization variables
nElemBeam1 = 20
nElemBeam2 = 10

# Range of angular velocity [rad/s]
ωRange = 2*π/60*[0,500,750]

# Range of beam tip angles [rad]
tipAngleRange = π/180*collect(0:2.5:45)

# Number of modes
nModes = 8

# Beams
stiffnessMatrix = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam1 = create_Beam(name="beam1",length=L1,nElements=nElemBeam1,C=[stiffnessMatrix],I=[inertiaMatrix])
beam2 = create_Beam(name="beam2",length=L2,nElements=nElemBeam2,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321")

# BCs
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Create model
sweptTipRotor = create_Model(name="sweptTipRotor",beams=[beam1,beam2],BCs=[clamp],initialPosition=r0,units=create_UnitsSystem(length="in",force="lbf",frequency="Hz"))

# Initialize outputs
numFreqs = Matrix{Vector{Float64}}(undef,length(ωRange),length(tipAngleRange))

# Loop over sweep variables
for (i,ω) in enumerate(ωRange)
    # Set angular velocity of basis A 
    sweptTipRotor.ω_A = [0; 0; ω]
    # Loop tip angles
    for (j,tipAngle) in enumerate(tipAngleRange)
        # Display progress
        ωRPM = round(Int,ω/(2*π/60))
        tipAngleDeg = round(tipAngle*180/π,digits=1)
        display("Solving for ω = $ωRPM rpm and tip angle = $tipAngleDeg deg")
        # Update beam2 angle with tip angle
        beam2.p0[1] = -tipAngle
        update_beam!(beam2)
        # Update model
        sweptTipRotor.beams = [beam1,beam2]
        update_model!(sweptTipRotor)
        # Create and solve eigenproblem
        global problem = create_EigenProblem(model=sweptTipRotor,nModes=nModes)
        solve!(problem)
        # Get outputs
        numFreqs[i,j] = problem.frequenciesOscillatory
    end
end

# Load experimental values
expTipAngles = [0, 15, 30, 45]
expFreqs1 = [1.4 1.8 1.7 1.6; 10.2 10.1 10.2 10.2; 14.8 14.4 14.9 14.7]
expFreqs2 = [10.3 10.2 10.4 10.4; 25.2 25.2 23.7 21.6; 36.1 34.8 30.7 26.1]
expFreqs3 = [27.7 27.2 26.6 24.8; 47.0 44.4 39.3 35.1; 62.9 55.9 48.6 44.8]
expFreqs4 = [95.4 87.5 83.7 78.8; 106.6 120.1 122.6 117.7; 132.7 147.3 166.2 162.0]

println("Finished sweptTipRotor.jl")
