# # Static structural analysis of the Pazy wing
# This example illustrates how to set up a static structural analysis, using the Technion's Pazy wing benchmark. More specifically, we'll simulate a coupled bending-torsion test. The sectional properties of the wing's spar and the aerodynamic tip loss function were defined by [Riso and Cesnik](https://doi.org/10.2514/6.2022-2313). That data is publicly available at [https://github.com/UM-A2SRL/AePW3-LDWG](https://github.com/UM-A2SRL/AePW3-LDWG).

#md # ![](../assets/pazyTorsionTest.png)
#md # *Pazy wing coupled bending-torsion test* by [Avin et al.](https://doi.org/10.2514/1.J060621)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/staticStructural/PazyWingTorsionTest.jl).

# ### Problem setup
# Let's begin by setting the variables of our problem. In this example we will analyze the displacements and twist of the clamped wing caused by the attachment of a range of tip masses. Notice that we bring into scope some fixed geometrical and discretization properties of the wing's beam through the function [`geometrical_properties_Pazy`](@ref geometrical_properties_Pazy).
using AeroBeams

## Fixed geometrical properties
nElem,L,chord,normSparPos,_ = geometrical_properties_Pazy()

## Set tip mass range
mRange = collect(0:0.1:3)

## Initialize outputs
tip_twist = Array{Float64}(undef,length(mRange))
tip_OOP = Array{Float64}(undef,length(mRange))
problem = Array{SteadyProblem}(undef,length(mRange))
#md nothing #hide

# ### Solving the problem
# In the following loop, we create new model instances for each tip mass value, create and solve the steady problem, and then extract the outputs of interest. The tip mass is positioned in such a way to cause not only bending but also significant torsion of the wing. The model creation process is streamlined with the function [`create_Pazy`](@ref create_Pazy), taking the appropriate inputs. Notice that the twist angle is computed as the arcsine of the difference in "height" between the leading edge and trailing edge of the wingtip (throught the rotation of the local basis vector pointing in the x2-direction).
## Sweep tip mass
for (i,m) in enumerate(mRange)
    ## Create model with current tip mass 80 mm behind the trailing-edge
    PazyWingTorsionTest,_ = create_Pazy(tipMass=m,ηtipMass=[0;-(chord*(1-normSparPos)+0.08);0])
    ## Create and solve problem
    problem[i] = create_SteadyProblem(model=PazyWingTorsionTest)
    solve!(problem[i])
    ## Get twist and OOP displacement at midchord
    tipNodeFinalState = problem[i].loadstepHistory.states[end][nElem].nodal
    tip_twist[i] = twist_angle(tipNodeFinalState.p2_b)
    tip_OOP[i] = tipNodeFinalState.u2[3] - chord*(1/2-normSparPos)*sind(tip_twist[i])
end
#md nothing #hide

# ### Post-processing
# The post-processing begins by loading the reference data.
## Load reference data
u3Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/torsion_u3VsMass_Exp.txt")
u3UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/torsion_u3VsMass_UMNAST.txt")
θExp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/torsion_thetaVsMass_Exp.txt")
θUMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/torsion_thetaVsMass_UMNAST.txt")
#md nothing #hide

#md # We can now plot the out-of-plane (OOP) tip displacement and tip twist as functions of the tip mass. The following "experimental" results were taken from Figure 17 of the paper by [Avin et al.](https://doi.org/10.2514/1.J060621). The correlation with the reference data is good.
#md gr()
#md ENV["GKSwstype"] = "100" #hide

## Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
#md plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP displacement offset [% semispan]", xlims=[mRange[1],mRange[end]])
#md plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=2, label="AeroBeams")
#md plot!(u3UMNAST[1,:], u3UMNAST[2,:], c=:blue, ls=:dash, lw=2, label="UM/NAST")
#md scatter!(u3Exp[1,:], u3Exp[2,:], mc=:red, ms=4, msw=0, label="Experimental")
#md savefig("PazyWingTorsionTest_OOP.svg") #hide
#md nothing #hide

## Tip twist vs. tip mass
#md plt2 = plot(xlabel="Tip mass [kg]", ylabel="Tip twist [deg]", xlims=[mRange[1],mRange[end]])
#md plot!(mRange, tip_twist, c=:black, lw=2, label="AeroBeams")
#md plot!(θUMNAST[1,:], θUMNAST[2,:], c=:blue, ls=:dash, lw=2, label="UM/NAST")
#md scatter!(θExp[1,:], θExp[2,:], mc=:red, ms=4, msw=0, label="Experimental")
#md savefig("PazyWingTorsionTest_twist.svg") #hide
#md nothing #hide

#md # ![](PazyWingTorsionTest_OOP.svg)
#md # ![](PazyWingTorsionTest_twist.svg)

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingTorsionTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem[end],view=(30,30),save=true,savePath=string(relPath,"/PazyWingTorsionTest_deformation.pdf"))
display(deformationPlot)

# Plot configurations
ts = 10
fs = 16
lw = 2
ms = 4
msw = 0
gr()

# Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,3], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=lw, label="AeroBeams")
plot!(u3UMNAST[1,:], u3UMNAST[2,:], c=:black, ls=:dot, lw=lw, label="UM/NAST - Riso & Cesnik (2022)")
scatter!(u3Exp[1,:], u3Exp[2,:], mc=:black, ms=ms, msw=msw, label="Exp. - Avin et al. (2022)")
display(plt1)
savefig(string(absPath,"/PazyWingTorsionTest_OOP.pdf"))

# Tip twist vs. tip mass
plt2 = plot(xlabel="Tip mass [kg]", ylabel="Tip twist [deg]", xlims=[0,3], tickfont=font(ts), guidefont=font(fs))
plot!(mRange, tip_twist, c=:black, lw=lw, label=false)
plot!(θUMNAST[1,:], θUMNAST[2,:], c=:black, ls=:dot, lw=lw, label=false)
scatter!(θExp[1,:], θExp[2,:], mc=:black, ms=ms, msw=msw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingTorsionTest_twist.pdf"))

println("Finished PazyWingTorsionTest.jl") #src