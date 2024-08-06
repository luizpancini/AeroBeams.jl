using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Pazy wing
wing,L,nElem,chord,normSparPos,_ = create_Pazy()

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingBendingTest = create_Model(name="PazyWingBendingTest",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.807])

# Set tip mass range, and initialize outputs
mRange = collect(0:0.1:3)
tip_OOP = Array{Float64}(undef,length(mRange))

# Sweep tip mass
for (i,m) in enumerate(mRange)
    # Display progress
    println("Solving for tip mass = $m kg")
    # Reset point inertias on beam
    wing.pointInertias = Vector{PointInertia}()
    update_beam!(wing)
    # Create tip mass (at midchord)
    tipMass = PointInertia(elementID=nElem,η=[L/nElem/2;-chord*(1/2-normSparPos);0],mass=m)
    # Add tip mass to beam and update model
    add_point_inertias_to_beam!(wing,inertias=[tipMass])
    update_model!(PazyWingBendingTest)
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
bending_u3VsMass_Exp = readdlm(string(pwd(),"/test/referenceData/Pazy/bending_u3VsMass_Exp.txt"))
bending_u3VsMass_UMNAST = readdlm(string(pwd(),"/test/referenceData/Pazy/bending_u3VsMass_UMNAST.txt"))

# Plots
# ------------------------------------------------------------------------------
# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath="/test/outputs/figures/PazyWingBendingTest/PazyWingBendingTest_deformation.pdf")
display(deformationPlot)
# Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP displacement offset [% semispan]", xlims=[0,3])
plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=2, label="AeroBeams")
plot!(bending_u3VsMass_UMNAST[1,:], bending_u3VsMass_UMNAST[2,:], c=:blue, ls=:dash, lw=2, label="UM/NAST")
scatter!(bending_u3VsMass_Exp[1,:], bending_u3VsMass_Exp[2,:], mc=:red, ms=3,msw=0, label="Exp.")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyWingBendingTest/PazyWingBendingTest_OOP.pdf"))

println("Finished PazyWingBendingTest.jl")