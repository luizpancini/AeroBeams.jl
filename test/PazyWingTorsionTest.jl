using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Pazy wing
wing,L,nElem,chord,normSparPos,_ = create_Pazy()

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingTorsionTest = create_Model(name="PazyWingTorsionTest",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.807])

# Set tip mass range, and initialize outputs
mRange = collect(0:0.1:3)
tip_twist = Array{Float64}(undef,length(mRange))
tip_OOP = Array{Float64}(undef,length(mRange))

# Sweep tip mass
for (i,m) in enumerate(mRange)
    # Display progress
    println("Solving for tip mass = $m kg")
    # Reset point inertias on beam
    wing.pointInertias = Vector{PointInertia}()
    update_beam!(wing)
    # Create tip mass (80 mm behind trailing-edge)
    tipMass = PointInertia(elementID=nElem,η=[L/nElem/2;-chord*(1-normSparPos)-0.08;0],mass=m)
    # Add tip mass to beam and update model
    add_point_inertias_to_beam!(wing,inertias=[tipMass])
    update_model!(PazyWingTorsionTest)
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
torsion_u3VsMass_Exp = readdlm(string(pwd(),"/test/referenceData/Pazy/torsion_u3VsMass_Exp.txt"))
torsion_u3VsMass_UMNAST = readdlm(string(pwd(),"/test/referenceData/Pazy/torsion_u3VsMass_UMNAST.txt"))
torsion_thetaVsMass_Exp = readdlm(string(pwd(),"/test/referenceData/Pazy/torsion_thetaVsMass_Exp.txt"))
torsion_thetaVsMass_UMNAST = readdlm(string(pwd(),"/test/referenceData/Pazy/torsion_thetaVsMass_UMNAST.txt"))

# Plots
# ------------------------------------------------------------------------------
# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath="/test/outputs/figures/PazyWingTorsionTest/PazyWingTorsionTest_deformation.pdf")
display(deformationPlot)
# Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
gr()
plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP displacement offset [% semispan]", xlims=[0,3])
plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=2, label="AeroBeams")
plot!(torsion_u3VsMass_UMNAST[1,:], torsion_u3VsMass_UMNAST[2,:], c=:blue, ls=:dash, lw=2, label="UM/NAST")
scatter!(torsion_u3VsMass_Exp[1,:], torsion_u3VsMass_Exp[2,:], mc=:red, ms=3,msw=0, label="Exp.")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyWingTorsionTest/PazyWingTorsionTest_OOP.pdf"))
# Tip twist vs. tip mass
plt2 = plot(xlabel="Tip mass [kg]", ylabel="Tip twist [deg]", xlims=[0,3])
plot!(mRange, tip_twist, c=:black, lw=2, label="AeroBeams")
plot!(torsion_thetaVsMass_UMNAST[1,:], torsion_thetaVsMass_UMNAST[2,:], c=:blue, ls=:dash, lw=2, label="UM/NAST")
scatter!(torsion_thetaVsMass_Exp[1,:], torsion_thetaVsMass_Exp[2,:], mc=:red, ms=3,msw=0, label="Exp.")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyWingTorsionTest/PazyWingTorsionTest_twist.pdf"))

println("Finished PazyWingTorsionTest.jl")