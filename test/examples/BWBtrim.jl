using AeroBeams, DelimitedFiles

# Stiffness factor
λ = 1e0

# Aerodynamic solver
aeroSolver = Indicial()

# Tip loss option
hasTipCorrection = true

# Model 
BWB = create_BWB(aeroSolver=aeroSolver,stiffnessFactor=λ,δElevIsTrimVariable=true,thrustIsTrimVariable=true,hasTipCorrection=hasTipCorrection)

# Set NR system solver 
relaxFactor = 0.5
displayStatus = false
maxiter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxiter,displayStatus=displayStatus)

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[3],nodesSides=[2],ku=ku,kp=kp)
add_springs_to_beam!(beam=BWB.beams[2],springs=[spring1])
add_springs_to_beam!(beam=BWB.beams[3],springs=[spring2])

# Set airspeed range and initialize outputs
URange = collect(30:5:160)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Trimming for airspeed U = $U m/s")
    # Update airspeed on model
    set_motion_basis_A!(model=BWB,v_A=[0;U;0])
    # Set initial guess solution as previous known solution
    x0 = i == 1 ? zeros(0) : problem.x
    # Create and solve trim problem
    global problem = create_TrimProblem(model=BWB,systemSolver=NR,x0=x0)
    solve!(problem)
    # Trim results
    trimAoA[i] = problem.aeroVariablesOverσ[end][BWB.beams[3].elementRange[1]].flowAnglesAndRates.αₑ*180/π
    trimThrust[i] = problem.x[end-1]*problem.model.forceScaling 
    trimδ[i] = problem.x[end]*180/π 
    println("AoA = $(trimAoA[i]), T = $(trimThrust[i]), δ = $(trimδ[i])")
end

# Load reference solution
trimAoARef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "BWB", "trimAoA.txt"))
trimThrustRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "BWB", "trimThrust.txt"))
trimδRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "BWB", "trimDelta.txt"))

println("Finished BWBtrim.jl")