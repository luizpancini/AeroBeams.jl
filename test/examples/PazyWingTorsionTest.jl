using AeroBeams, DelimitedFiles

# Fixed geometrical properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Set tip mass range, and initialize outputs
mRange = collect(0:0.1:3)
tip_twist = Array{Float64}(undef,length(mRange))
tip_OOP = Array{Float64}(undef,length(mRange))

# Sweep tip mass
for (i,m) in enumerate(mRange)
    # Create model with current tip mass 80 mm behind the trailing-edge
    PazyWingTorsionTest,_ = create_Pazy(tipMass=m,ξtipMass=[0;-(chord*(1-normSparPos)+0.08);0])
    # Create and solve problem
    global problem = create_SteadyProblem(model=PazyWingTorsionTest)
    solve!(problem)
    # Get OOP displacement at midchord
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist[i] = asind(Δ[3])
    tip_OOP[i] = problem.nodalStatesOverσ[end][nElem].u_n2[3] - chord*(1/2-normSparPos)*sind(tip_twist[i])
end

# Load reference data
u3Exp = readdlm("test/referenceData/Pazy/torsion_u3VsMass_Exp.txt")
u3UMNAST = readdlm("test/referenceData/Pazy/torsion_u3VsMass_UMNAST.txt")
θExp = readdlm("test/referenceData/Pazy/torsion_thetaVsMass_Exp.txt")
θUMNAST = readdlm("test/referenceData/Pazy/torsion_thetaVsMass_UMNAST.txt")

println("Finished PazyWingTorsionTest.jl")