using AeroBeams

# Fixed geometrical properties
nElem,L,chord,normSparPos,_ = geometrical_properties_Pazy()

# Set tip mass range, and initialize outputs
mRange = collect(0:0.1:3)
tip_OOP = Vector{Float64}(undef,length(mRange))
problem = Vector{SteadyProblem}(undef,length(mRange))

# Sweep tip mass
for (i,m) in enumerate(mRange)
    # Display progress
    println("Solving for tip mass = $m kg")
    # Create model with current tip mass at the trailing-edge
    PazyWingBendingTest,_ = create_Pazy(aeroSolver=QuasiSteady(),tipMass=m,ηtipMass=[0; -chord*(1-normSparPos); 0])
    # Create and solve problem
    problem[i] = create_SteadyProblem(model=PazyWingBendingTest,systemSolver=create_NLsolver())
    solve!(problem[i])
    # Get OOP displacement at midchord
    tipNodeFinalState = problem[i].loadstepHistory.states[end][nElem].nodal
    tip_twist = twist_angle(tipNodeFinalState.p2_b)
    tip_OOP[i] = tipNodeFinalState.u2[3] - chord*(1/2-normSparPos)*sind(tip_twist)
end

# Load reference data
u3Exp = readdlm(pkgdir(AeroBeams,"test","referenceData","Pazy", "bending_u3VsMass_Exp.txt"))
u3UMNAST = readdlm(pkgdir(AeroBeams,"test","referenceData","Pazy", "bending_u3VsMass_UMNAST.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingBendingTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape at the final tip mass
deformationPlot = plot_steady_deformation(problem[end],view=(30,30),save=true,savePath=string(relPath,"/PazyWingBendingTest_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
gr()

# Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP displacement offset [% semispan]", xlims=[0,3])
plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=lw, label="AeroBeams")
plot!(u3UMNAST[1,:], u3UMNAST[2,:], c=:blue, ls=:dash, lw=lw, label="UM/NAST")
scatter!(u3Exp[1,:], u3Exp[2,:], mc=:red, ms=ms, msw=msw, label="Exp.")
display(plt1)
savefig(string(absPath,"/PazyWingBendingTest_OOP.pdf"))

println("Finished PazyWingBendingTest.jl")