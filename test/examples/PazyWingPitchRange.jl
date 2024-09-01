using AeroBeams, DelimitedFiles

# Aerodynamic solver
aeroSolver = Indicial()

# Derivation method
derivationMethod = AD()

# Flag for upright position
upright = true

# Set root angle and airspeed ranges, and initialize outputs
θRange = [3;5;7]
URange = collect(0:1:60)
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
tip_IP = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))
tip_AoA = Array{Float64}(undef,length(θRange),length(URange))

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Update model
        PazyWingPitchRange,nElem,L,chord,normSparPos = create_Pazy(aeroSolver=aeroSolver,derivationMethod=derivationMethod,upright=upright,θ=θ*π/180,airspeed=U)
        # Create and solve problem
        global problem = create_SteadyProblem(model=PazyWingPitchRange)
        solve!(problem)
        # Get tip twist, AoA and OOP displacement at midchord
        tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist[i,j] = asind(Δ[3])
        tip_OOP[i,j] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist[i,j]))
        tip_IP[i,j] = -problem.nodalStatesOverσ[end][nElem].u_n2[2]
        tip_AoA[i,j] = problem.model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
    end
end

# Load reference data
tip_u3VsU_rootPitch5 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "Pazy", "tip_u3VsU_rootPitch5.txt"))
tip_u3VsU_rootPitch7 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "Pazy", "tip_u3VsU_rootPitch7.txt"))

println("Finished PazyWingPitchRange.jl")