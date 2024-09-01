using AeroBeams, DelimitedFiles

# Fixed geometrical properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Set tip mass range, and initialize outputs
mRange = collect(0:0.1:3)
tip_OOP = Array{Float64}(undef,length(mRange))

# Sweep tip mass
for (i,m) in enumerate(mRange)
    # Display progress
    println("Solving for tip mass = $m kg")
    # Create model with current tip mass at the trailing-edge
    PazyWingBendingTest,_ = create_Pazy(tipMass=m,ξtipMass=[0;-chord*(1-normSparPos);0])
    # Create and solve problem
    global problem = create_SteadyProblem(model=PazyWingBendingTest)
    solve!(problem)
    # Get OOP displacement at midchord
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist = asind(Δ[3])
    tip_OOP[i] = problem.nodalStatesOverσ[end][nElem].u_n2[3] - chord*(1/2-normSparPos)*sind(tip_twist)
end

# Load reference data
u3Exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "Pazy", "bending_u3VsMass_Exp.txt"))
u3UMNAST = readdlm(joinpath(dirname(@__DIR__), "referenceData", "Pazy", "bending_u3VsMass_UMNAST.txt"))

println("Finished PazyWingBendingTest.jl")