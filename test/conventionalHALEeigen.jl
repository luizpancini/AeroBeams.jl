using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Stiffness factor 
λ = 10

# Model
conventionalHALE,leftWing,_ = create_conventional_HALE(aeroSolver=Inflow(),derivationMethod=AD(),nElemWing=20,nElemTailBoom=10,nElemHorzStabilizer=10,stiffnessFactor=λ,stabilizersAero=true,includeVS=false,wingCd0=1e-2,k2=0)

# Trim thrust force 
trimThrust = create_BC(name="thrust",beam=leftWing,node=leftWing.nodeRange[end],types=["Ff2A"],values=[0],toBeTrimmed=trues(1))
conventionalHALE.BCs = [trimThrust]

# Set airspeed on model
U = 20
set_motion_basis_A!(model=conventionalHALE,v_A=[0;U;0])

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(initialLoadFactor=0.5,maximumIterations=100)

# Solve trim problem and get solution
trimproblem = create_TrimProblem(model=conventionalHALE,systemSolver=NR)
solve!(trimproblem)
x0 = trimproblem.x[1:end-1]
trimThrust = trimproblem.x[end] * conventionalHALE.forceScaling

# Set thrust force and update model
thrust = create_BC(name="thrust",beam=leftWing,node=leftWing.nodeRange[end],types=["Ff2A"],values=[trimThrust])
conventionalHALE.BCs = [thrust]
conventionalHALE.skipValidationMotionBasisA = true
update_model!(conventionalHALE)

# Number of modes and frequency bounds
nModes = 12
frequencyLimits = [1e-2, 200]

# Create and solve eigenproblem
eigenproblem = create_EigenProblem(model=conventionalHALE,nModes=nModes,frequencyFilterLimits=frequencyLimits,x0=x0)
solve!(eigenproblem)

# Get frequency and damping of rigid body modes
freqs = eigenproblem.frequenciesOscillatory
damps = eigenproblem.dampingsOscillatory

for (i,freq,damp) in zip(1:length(freqs),freqs,damps)
    println("Eigenvalue $i: $damp +/- $freq i")
end

# # Compare with reference solution
# ωₚRef, ωₛₚRef = 0.319, 5.67
# ζₚRef, ζₛₚRef = 0.0709, 0.905

# ϵωₚ = freqs[1]/ωₚRef - 1
# ϵωₛₚ = freqs[2]/ωₛₚRef - 1
# ϵζₚ = damps[1]/ζₚRef - 1
# ϵζₛₚ = damps[2]/ζₛₚRef - 1

# println("Relative errors")
# println("Phugoid frequency: $ϵωₚ")
# println("Short period frequency: $ϵωₛₚ")
# println("Phugoid damping: $ϵζₚ")
# println("Short period damping: $ϵζₛₚ")

println("Finished conventionalHALEeigen.jl")