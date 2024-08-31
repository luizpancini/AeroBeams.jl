using AeroBeams, DelimitedFiles

# Aerodynamic solver
aeroSolver = Indicial()

# Option for stabilizers
stabilizersAero = false
includeVS = false

# Parasite drag coefficients
wingCd0 = stabsCd0 = 0

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10

# Wing precurvature
k2 = 0.0

# NR system solver 
maxit = 100
displayStatus = true
NR = create_NewtonRaphson(maximumIterations=maxit,displayStatus=displayStatus)

# Set stiffness factor and airspeed ranges, and initialize outputs
λRange = [1; 50]
URange = collect(20:1:35)
trimAoA = Array{Float64}(undef,length(λRange),length(URange))
trim_u1 = Array{Vector{Float64}}(undef,length(λRange),length(URange))
trim_u3 = Array{Vector{Float64}}(undef,length(λRange),length(URange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Model
    conventionalHALE,leftWing,rightWing,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,k2=k2)
    # Get element ranges and nodal arclength positions of right wing
    global x1_0 = vcat([vcat(rightWing.elements[e].r_n1[1],rightWing.elements[e].r_n2[1]) for e in 1:rightWing.nElements]...)
    global x3_0 = vcat([vcat(rightWing.elements[e].r_n1[3],rightWing.elements[e].r_n2[3]) for e in 1:rightWing.nElements]...)
    rWGlobalElemRange = rightWing.elementRange[1]:rightWing.elementRange[end]
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Trimming for λ = $λ, U = $U m/s")
        # Update airspeed on model
        set_motion_basis_A!(model=conventionalHALE,v_A=[0;U;0])
        # Set initial guess solution as previous known solution
        x0 = j == 1 ? zeros(0) : problem.x
        # Create and solve trim problem
        global problem = create_TrimProblem(model=conventionalHALE,systemSolver=NR,x0=x0)
        solve!(problem)
        # Trim results
        trimAoA[i,j] = problem.aeroVariablesOverσ[end][rWGlobalElemRange[1]].flowAnglesAndRates.αₑ*180/π
        trim_u1[i,j] = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in rWGlobalElemRange]...)
        trim_u3[i,j] = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in rWGlobalElemRange]...)
        println("AoA = $(trimAoA[i,j]) deg")
    end
end

# Load reference solutions
trimAoAERef = readdlm("test/referenceData/conventionalHALE/trimAoAVsAirspeedElastic.txt")
trimAoARRef = readdlm("test/referenceData/conventionalHALE/trimAoAVsAirspeedRigid.txt")
trimDispRef = readdlm("test/referenceData/conventionalHALE/trimDispAtU25.txt")

println("Finished conventionalHALEtrim.jl")