using AeroBeams, DelimitedFiles

# Wing airfoil
# wingAirfoil = NACA23012A
wingAirfoil = HeliosWingAirfoil

# Option for reduced chord
reducedChord = false

# TF to include beam pods 
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Set NR system solver 
relaxFactor = 0.5
displayStatus = false
NR = create_NewtonRaphson(ρ=relaxFactor,displayStatus=displayStatus)

# Airspeed
U = 40*0.3048

# Set stiffness factor and payload ranges, and initialize outputs
λRange = [1,50]
PRange = collect(0:20:500)
problem = Array{TrimProblem}(undef,length(λRange),length(PRange))
trimAoA = Array{Float64}(undef,length(λRange),length(PRange))
trimThrust = Array{Float64}(undef,length(λRange),length(PRange))
trimδ = Array{Float64}(undef,length(λRange),length(PRange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Trimming for λ = $λ payload = $P lb")
        # Model and its beams
        helios,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,reducedChord=reducedChord,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        # Set initial guess solution as previous known solution
        x0 = (j==1) ? zeros(0) : problem[i,j-1].x
        # Create and solve trim problem[i,j]
        problem[i,j] = create_TrimProblem(model=helios,systemSolver=NR,x0=x0)
        solve!(problem[i,j])
        # Trim results
        trimAoA[i,j] = problem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
        trimThrust[i,j] = problem[i,j].x[end-1]*problem[i,j].model.forceScaling
        trimδ[i,j] = problem[i,j].x[end]*180/π
        println("AoA = $(trimAoA[i,j]), T = $(trimThrust[i,j]), δ = $(trimδ[i,j])")
    end
end

# Load reference data
αFlexibleRef = readdlm("test/referenceData/Helios/trim_AoA_flexible.txt")
αRigidRef = readdlm("test/referenceData/Helios/trim_AoA_rigid.txt")
δFlexibleRef = readdlm("test/referenceData/Helios/trim_delta_flexible.txt")
δRigidRef = readdlm("test/referenceData/Helios/trim_delta_rigid.txt")

println("Finished heliosTrim.jl")