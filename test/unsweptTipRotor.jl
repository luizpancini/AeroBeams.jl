using AeroBeams, LinearAlgebra, Plots, ColorSchemes

## User inputs (problem definition)
#-------------------------------------------------------------------------------
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

## Problem setup
#-------------------------------------------------------------------------------
# Beams
complianceMatrix = zeros(6,6)
complianceMatrix[1:3,1:3] .= [   1/(E11*A) -ν12/(G12*A*Ksy)  -ν13/(G13*A*Ksz);
                          -ν12/(G12*A*Ksy)    1/(G12*A*Ksy)      -ν23/(G23*A);
                          -ν13/(G13*A*Ksz)     -ν23/(G23*A)     1/(G13*A*Ksz)]
complianceMatrix[4:6,4:6] .= diagm([1/(G23*J*Kt); 1/(E11*Iy); 1/(E11*Iz)])
stiffnessMatrix = inv(complianceMatrix)
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Create model
unsweptTipRotor = create_Model(name="unsweptTipRotor",beams=[beam],BCs=[clamp],initialPosition=r0)

# Initialize outputs
numFreqs = Vector{Vector{Float64}}(undef,length(ωRange))

# Loop over angular velocities
for (i,ω) in enumerate(ωRange)
    # Set angular velocity of basis A and update model
    unsweptTipRotor.ω_A = [0; 0; ω]
    update_model!(unsweptTipRotor)
    # Display progress
    display("Solving for ω = $(ωRPM[i]) rpm")
    # Create and solve eigenproblem
    global problem = create_EigenProblem(model=unsweptTipRotor,nModes=nModes)
    solve!(problem)
    # Get outputs
    numFreqs[i] = problem.frequenciesOscillatory
end

## Plots
#-------------------------------------------------------------------------------
# Load experimental values from Epps & Chandra (1996)
ωRPMExp = [0,250,500,750]
expFreqs = [5.2374e+00   7.4239e+00   1.1018e+01   1.5820e+01;
            3.2999e+01   3.4582e+01   4.0188e+01   4.7001e+01;
            9.0439e+01   9.1112e+01   9.7733e+01   1.0558e+02;
            1.7237e+02   1.7532e+02   1.7967e+02   1.8874e+02]

# Set colormap and labels
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, 4))
labels = ["1B", "2B", "3B", "4B"]

# Plot frequency over angular velocity for several modes
plt1 = plot(xlabel="Angular velocity [rpm]", ylabel="Frequency [Hz]", title="Bending modes", xticks=ωRPMExp, yticks=collect(0:20:200))
modes = [1,2,5,6]
plot!([NaN], [NaN], lc=:black, lw=2, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=5, label="Epps & Chandra (1996)")
for (m,mode) in enumerate(modes)  
    numFreqsMode = [numFreqs[j][mode]/(2*π) for j in 1:length(ωRPM)]
    plot!(ωRPM,numFreqsMode, lc=colors[m], lw=2, label=false)
    scatter!(ωRPMExp,expFreqs[m,:], mc=colors[m], ms=5, label=false)
    plot!([NaN], [NaN], lc=colors[m], m=colors[m], lw=2, ms=5, label=labels[m])
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/unsweptTipRotor.pdf"))

println("Finished unsweptTipRotor.jl")